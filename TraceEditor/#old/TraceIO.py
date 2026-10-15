import configparser
import datetime
import os
import re
import time
import openpyxl
from pathlib import Path

import h5py
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import uncertainties as unc
from lmfit import Model
from mpl_toolkits.axes_grid1 import Divider, Size
from nptdms import TdmsFile
import Profiler as pr
from tqdm import tqdm

import ForceSpectroscopy as fs

"""
This Module is to simplify handling data files and to enhance synergy between different projects.
The basis is that all data consists of (time-)traces. Such a trace is defined by the variables:
filename, channel, label. For example:
in FORCE SPECTROSCOPY a label refers to a particular bead, a channel would be the extension, 'Z (um)'.
in SINGLE-MOLECULE FLUORESCENCE MICROSCOPY a label refers to a particular molecule, a channel would be it's
intensity 'I (kHz)', or its position 'X (um)'.
in FLUORESCENCE CORRELATION SPECTROSCOPY a label refers to correlation signal, i.e. 'R632xG524'. The channel would be
it's correlation curve, 'G-1 (-)'.
Shared traces ('Time (s)', 'F (pN)', 'Tau (s)', ...) do not have a label (or to be more precise, they have the
label ' shared', which is the default label and can be left out).

All traces are stored in the same HDF5 file, and are organized in groups. HDF5 files (extension *.h5) are
hierarchical binary files, whose structure has been amply described (https://www.hdfgroup.org/). Download and
install HDFview (https://www.hdfgroup.org/downloads/hdfview/) to directly see the file contents.

Access traces in Python by h5_read_trace() and h5_write_trace(). It's good practice to maintain the following file
structure: ../users/USERNAME/data/DATE/FILE.h5
A trace is then uniquely defined by: ../users/USERNAME/data/YYYYMMDDD/FILE.h5 /traces/label > channel

In addition to (time-)traces there may be one ore more parameters that are linked to a particular label, like offsets, 
start times, fit parameters, mean values and.or standard deviations of a trace, etc. These parameters are also
stored in the HDF5 file and can be read/written with the functions h5_read_par()/h5_write_par().

The function h5_contents() lists all channels, labels and parameters in the file.

Fitted curves are NOT stored in the file. However, they can be retrieved from the stored parameters using
h5_read_model(). You can constrain one ore more parameters to values that are not in the HDF5 file using constraints.
When not given in the file or as constraint, the function will use default values, as defined in the fit-function.
You can easily add new fit functions as long as you keep the standerd format of a fit function (look up an example 
like 'ExponentialDrift'). 

A dataframe of all stored parameters of the file can be obtained with the function h5_read_par_table().
Familiarize yourself with pandas (https://pandas.pydata.org/pandas-docs/stable/getting_started/10min.html),
to achieve the most efficient downstream processing. It has similar functionality as Excel, and much more...

"""


class H5_data(object):
    def __init__(self, filename=None, label=None):
        self.contents = pd.DataFrame()

        self.filename = None
        self.label = None
        self.parameters = pd.DataFrame()
        self.traces = pd.DataFrame()
        self.shared_traces = []
        if filename is not None and label is not None:
            self.open(filename, label)

    def open(self, filename, label, save=True):
        if filename != self.filename or label != self.label:
            if save:
                self.save()
            if filename and Path(filename).exists():
                self.filename = filename
                self.label = label
                self.traces = pd.DataFrame()
                self.parameters = pd.DataFrame(columns=["value", "error"])
                self.shared_traces = []
                with h5py.File(filename, "r") as f:
                    for name, item in f[f"/traces/ shared"].items():
                        self.traces[name] = item[()]
                        self.shared_traces.append(name)
                    for name, item in f[f"/traces/{label}"].items():
                        self.traces[name] = item[()]
                    for name, item in f[f"/parameters/values"].items():
                        self.parameters.loc[name, "value"] = item[()][int(label)]
                    for name, item in f[f"/parameters/errors"].items():
                        self.parameters.loc[name, "error"] = item[()][int(label)]

                    labels = f["/parameters/label"][()][:-1].astype("int")
                    selected = f["/parameters/values/ selected"][()].astype("int")

                # self.parameters = pd.DataFrame(columns=['value', 'error']) ###>>>>>> to be removed
                df_contents = pd.DataFrame(
                    {
                        "filenames": [filename] * len(labels),
                        "labels": labels,
                        "selected": selected,
                    }
                )
                self.contents = pd.concat(
                    [self.contents, df_contents], ignore_index=True
                ).drop_duplicates()
                self.contents = self.contents.sort_values(by=["filenames", "labels"])
            else:
                print(f"@H5_data.open: File does not exist: {filename}")
                self.filename = None
                self.label = None

    def close(self):
        self.contents = self.contents[self.contents["filenames"] != self.filename]
        try:
            self.open(
                self.contents["filenames"].iloc[-1], self.contents["labels"].iloc[-1]
            )
        except IndexError:
            self.save()
            self.filename = None
            self.label = None

    def save(self, params_only=True):
        if not self.filename or not self.label:
            # print("Filename or label not set. Exiting.")
            return

        # stack = inspect.stack()
        # print(
        #     f'>>>>>>>Saving {self.filename}...{self.label}, Called by {stack[2].function}')
        if Path(filename).suffix == ".h5":
            with h5py.File(self.filename, "a") as f:
                # Save traces
                if not params_only:
                    for col in self.traces.columns:
                        shared_group = f.require_group(f"/traces/shared")
                        trace_group = f.require_group(f"/traces/{self.label}")
                        if col in shared_group:
                            del shared_group[col]  # Delete existing dataset
                        if col in trace_group:
                            del trace_group[col]  # Delete existing dataset
                        if col in self.shared_traces:
                            shared_group.create_dataset(
                                col, data=self.traces[col].to_numpy()
                            )
                        else:
                            trace_group.create_dataset(
                                col, data=self.traces[col].to_numpy()
                            )

                # Save parameters
                for col in ["values", "errors"]:
                    param_group = f.require_group(f"/parameters/{col}")
                    for idx, row in self.parameters.iterrows():
                        if idx in param_group:
                            # Overwrite only the corresponding label index if dataset exists
                            param_group[idx][int(self.label)] = row[col[:-1]]
                        else:
                            # Create a new dataset and fill with NaN or appropriate value for new label
                            # assuming label starts from 0
                            n_labels = int(self.label) + 1
                            new_data = np.full((n_labels,), np.nan)
                            new_data[int(self.label)] = row[col[:-1]]
                            param_group.create_dataset(idx, data=new_data)

    def read_xlsx(self, label, filename=None):
        if filename is None:
            filename = self.filename
        filename = Path(self.filename).with_suffix(".xlsx")

        traces = pd.read_excel(filename, sheet_name="shared", index_col=0)
        self.shared_traces = traces.columns
        self.traces = pd.read_excel(
            filename, sheet_name=str(label), index_col=0
        ).combine_first(traces)

        self.parameters = pd.DataFrame(columns=["value", "error"])
        df = pd.read_excel(filename, sheet_name="parameters", index_col=0).T[self.label]
        value_series = df[~df.index.str.startswith("e:")]
        error_series = df[df.index.str.startswith("e:")].rename(lambda x: x[2:])
        self.parameters = pd.DataFrame({"value": value_series, "error": error_series})

        self.filename = filename
        self.label = label
        return

    def save_xlsx(self):
        def strip_prefix(col):
            parts = col.split(":")
            if len(parts) > 1:
                return parts[-1]
            return col

        filename = Path(self.filename).with_suffix(".xlsx")

        pars_out = self.parameters[["value"]].rename({"value": self.label}, axis=1).T
        errs_out = (
            self.parameters[["error"]].dropna().rename({"error": self.label}, axis=1).T
        )
        errs_out.rename(columns=lambda x: "e:" + x, inplace=True)
        pars_out = pd.concat([pars_out, errs_out], axis=1)

        if not os.path.exists(filename):
            with pd.ExcelWriter(filename, engine="openpyxl") as writer:
                sorted_cols = sorted(pars_out.columns, key=strip_prefix)
                pars_out = pars_out[sorted_cols]
                pars_out.to_excel(writer, sheet_name="parameters")
        else:
            old_df = pd.read_excel(filename, sheet_name="parameters", index_col=0)
            pars_out = pars_out.combine_first(old_df)
            sorted_cols = sorted(pars_out.columns, key=strip_prefix)
            pars_out = pars_out[sorted_cols]
            with pd.ExcelWriter(
                filename, engine="openpyxl", mode="a", if_sheet_exists="replace"
            ) as writer:
                pars_out.to_excel(writer, sheet_name="parameters")

        with pd.ExcelWriter(
            filename, engine="openpyxl", mode="a", if_sheet_exists="replace"
        ) as writer:
            for sheet in ["shared", "traces"]:
                if sheet == "traces":
                    sheet = str(self.label)
                    new_df = self.traces[
                        self.traces.columns[
                            ~self.traces.columns.isin(self.shared_traces)
                        ]
                    ]
                else:
                    new_df = self.traces[self.shared_traces]
                if sheet in writer.book.sheetnames:
                    old_df = pd.read_excel(filename, sheet_name=sheet, index_col=0)
                    new_df = new_df.combine_first(old_df)
                new_df.to_excel(writer, sheet_name=sheet)

    def correct_z(self, z_in, drift=True, bead_offset=True):
        if bead_offset:
            if "Zbead (um)" in self.parameters.index:
                if pd.notna(self.parameters.loc["Zbead (um)", "value"]):
                    z_in -= self.parameters.loc["Zbead (um)", "value"]
        if drift:
            z_in -= h5_read_model(self, "ExponentialDrift")[0]
        return z_in


DELIMITER = " > "


def h5_contents(filename, to_file=False, label=None):
    with h5py.File(filename, "a") as f:
        if label is None:
            labels = sorted([i[0] for i in list(f["traces"].items())], key=natural_keys)
        else:
            labels = [label]

        channels = []
        for l in labels:
            channels += [i[0] for i in list(f[f"traces/{l}"].items())]
        channels = sorted(set(channels), key=natural_keys)
        try:
            parameters = sorted(
                list(set([i[0] for i in list(f[r"parameters/values"].items())]))
            )
        except KeyError:
            parameters = []

        if to_file:
            entry = rf"/parameters/label"
            try:
                stored_data = f[entry]
                stored_data[...] = np.asarray(labels).astype("S")
            except KeyError:
                f.create_dataset(entry, data=np.asarray(labels).astype("S"))
    try:
        labels.remove(" shared")
    except:
        pass
    return channels, labels, parameters


def h5_write_trace(trace, data, channel):
    data.traces[channel] = trace
    return


def h5_read_trace(data, channel):
    try:
        return data.traces[channel].to_numpy()
    except KeyError:
        return None


def h5_read_model(data, modelname, constraints={}):
    model, parameters, names = prepare_fit(modelname)

    for key, value in parameters.valuesdict().items():
        if key in list(constraints.keys()):
            parameters[key].value = constraints[key]
        else:
            if names[key] in data.parameters.index:
                parameters[key].value = data.parameters.loc[names[key], "value"]

    x_data = h5_read_trace(data, names["x"])
    y_model = model.eval(parameters, x=x_data)
    return y_model, parameters


def h5_write_par(data, name, value, error=0, type="local"):
    data.parameters.loc[name] = [value, error]
    return


def h5_read_par(data, name):
    par = list(data.parameters.loc[name])
    types = ["fit", "local", "global"]
    try:
        return par[0], par[1], types[int(par[2])]
    except (ValueError, IndexError):
        return 0, 0, 2


def h5_read_par_table(filename, par_names=None, type="values", selected=True):
    if not os.path.exists(filename):
        print(f"@h5_read_par: File does not exist: {filename}")
        return

    with h5py.File(filename, "r") as hf:
        df = pd.DataFrame(hf[rf"/parameters/label"][:-1], columns=["label"])
        df["label"] = df["label"].str.decode("utf-8")

        if selected:
            selection = np.array(hf[rf"/parameters/values/ selected"]) == 1
            df = df[selection]
        else:
            selection = np.ones(len(df)) == 1

        if par_names is None:
            par_names = sorted(list(hf[rf"/parameters/{type}"].keys()))
        for par_name in par_names:
            df[par_name] = hf[rf"/parameters/{type}/{par_name}"][selection]
    df = df.set_index("label")
    return df


def h5_list_pars(data, selection=None):
    parameters = data.parameters.index.sort_values()
    if selection is not None:
        parameters = [p for p in selection if p in parameters]
    types = ["fit", "local", "global"]
    pars_list = []
    for par in parameters:
        value = data.parameters.loc[par, "value"]
        try:
            type = types[int(data.parameters.loc[par, "type"])]
        except:
            type = "local"
        error = data.parameters.loc[par, "error"]
        if type == "fit":
            pars_list.append([par, round_significance(value, error), type])
        else:
            pars_list.append([par, f"{value:.6g}", type])
    return pars_list


def prepare_fit(func):
    """
    Prepare all settings for fitting or evaluating a function.
    Make sure the that the function:

    - has a fully filled out doc string
    - the independent variable is called 'x' and comes first
    - each variable, including the independent variable and the return entry, has a proper, readable name.
    - the readable name includes units in brackets
    - the readable name is terminated by a hash and optionally supplemented by a description.
      (i.e. a valid entry would be:  ':param x: F (pN) # The force that was applied.')

    :param func: a function that is implemented in an imported module
    :return: an LMfit model, the corresponding LMfit parameters, and a dictionary of the names
             that describe the used function variables.
    """
    model = Model(getattr(fs, func))
    params = model.make_params()

    doc = getattr(fs, func).__doc__

    names = {}
    names["x"] = doc.split("x:")[1].split("#")[0].strip()
    names["y"] = doc.split("return:")[1].split("#")[0].strip()
    for key, item in params.items():
        names[key] = doc.split(f"{key}:")[1].split("#")[0].strip()
        params[key].vary = False
        try:
            range = (
                doc.split(f"{key}:")[-1]
                .split("\n")[0]
                .split("]")[0]
                .split("[")[-1]
                .strip()
                .split(", ")
            )
            range = np.asarray(range).astype(float)
            params[key].set(min=range[0], max=range[1])
        except ValueError:
            pass
    return model, params, names


def timeit(method):
    def timed(*args, **kw):
        ts = time.time()
        result = method(*args, **kw)
        te = time.time()
        if "log_time" in kw:
            name = kw.get("log_name", method.__name__.upper())
            kw["log_time"][name] = int((te - ts) * 1000)
        else:
            print(f"{time.asctime()} @timeit: {method.__name__} took {te - ts:.3f} s")
        return result

    return timed


def natural_keys(text):
    """
    alist.sort(key=natural_keys) sorts in human order
    http://nedbatchelder.com/blog/200712/human_sorting.html
    (See Toothy's implementation in the comments)
    float regex comes from https://stackoverflow.com/a/12643073/190597
    """

    def atof(text):
        try:
            retval = float(text)
        except ValueError:
            retval = text
        return retval

    return [atof(c) for c in re.split(r"[+-]?([0-9]+(?:[.][0-9]*)?|[.][0-9]+)", text)]


def sort_list(alist):
    return alist.sort(key=natural_keys)


def round_significance(x, stderr=None, par_name=None, delimiter=" ± "):
    """
    Compute proper significant digits. Standard error has 1 significant digit, unless
    that digit equals 1. In that cas it has 2 significant digits. Accuracy of x value
    matches accuracy of error.
    :param x: value
    :param stderr: standard error or standard deviation
    :param par_name: name of variable, preferably with units in brackets: 'x (a.u.)'
    :return: If par_name is provided a string will result: 'x (a.u.) = value +/- err'.
             Otherwise the rounded values of x and stderr.
    """
    if pd.isna(x):
        str_value = "-"
    elif pd.isna(stderr) or stderr <= 0:
        str_value = f"{x:g}"
    else:
        position = -int(np.floor(np.log10(abs(stderr))))
        first_digit = int(stderr * 10**position)
        if first_digit == 1:
            position += 1
        x_str = "{:.{}f}".format(x, position)
        stderr_str = "{:.{}f}".format(stderr, position)
        str_value = f"{x_str}{delimiter}{stderr_str}"

    if par_name is not None:
        str_value = f"{par_name} = " + str_value

    return str_value


def fix_axes(axew=7, axeh=7):
    # axew = axew / 2.54
    # axeh = axeh / 2.54

    topmargin = 0.5

    # lets use the tight layout function to get a good padding size for our axes labels.
    fig = plt.gcf()
    fig.set_size_inches(axew + 1.5, axeh + topmargin)

    ax = plt.gca()
    fig.tight_layout()
    # obtain the current ratio values for padding and fix size
    oldw, oldh = fig.get_size_inches()
    l = ax.figure.subplotpars.left
    r = ax.figure.subplotpars.right
    t = ax.figure.subplotpars.top + topmargin * 1.5
    b = ax.figure.subplotpars.bottom + topmargin * 1.5

    # work out what the new  ratio values for padding are, and the new fig size.
    neww = axew + oldw * (1 - r + l)
    newh = axeh + oldh * (1 - t + b)
    newr = r * oldw / neww
    newl = l * oldw / neww
    newt = t * oldh / newh
    newb = b * oldh / newh

    # right(top) padding, fixed axes size, left(bottom) pading
    hori = [Size.Scaled(newr), Size.Fixed(axew), Size.Scaled(newl)]
    vert = [Size.Scaled(newt), Size.Fixed(axeh), Size.Scaled(newb)]

    divider = Divider(fig, (0.0, 0.0, 0.95, 0.95), hori, vert, aspect=False)
    # the width and height of the rectangle is ignored.

    ax.set_axes_locator(divider.new_locator(nx=1, ny=1))

    # we need to resize the figure


def key_press_action(event):
    print(event.key)
    fig = plt.gcf()
    if event.key == "z":
        print(fig.texts[0].get_text())
    return


def format_plot(
    xtitle="x (a.u.)",
    ytitle="y (a.u.)",
    title="",
    xrange=None,
    yrange=None,
    ylog=False,
    xlog=False,
    scale_page=1.0,
    aspect=0.5,
    save=None,
    boxed=True,
    GUI=False,
    ref="",
    legend=None,
    fig=None,
    ax=None,
    txt=None,
    empty_canvas=False,
):
    # adjust the format to nice looks
    from matplotlib.ticker import AutoMinorLocator
    import os
    import subprocess

    page_width = 7  # inches ( = A4 width - 1 inch margin)
    margins = (0.55, 0.45, 0.2, 0.2)  # inches
    fontsize = 14

    if fig == None:
        fig = plt.gcf()
    if ax == None:
        ax = plt.gca()

    if empty_canvas:
        ax.set_axis_off()
        ax.set_xticks([])
        ax.set_yticks([])
        return fig, ax

    # Set up figure
    fig_width = page_width * scale_page * 0.95
    fig_height = (
        (fig_width - (margins[0] + margins[2])) * aspect + margins[1] + margins[3]
    )

    if txt:
        txt_width = 2
    else:
        txt_width = 0

    fig.set_size_inches(fig_width + txt_width, fig_height)

    # Set up axes
    ax_rect = [margins[0] / fig_width]
    ax_rect.append(margins[1] / fig_height)
    ax_rect.append((fig_width - margins[0] - margins[2]) / (fig_width + txt_width))
    ax_rect.append((fig_height - margins[1] - margins[3]) / fig_height)

    # ax_rect = [margins[0] / fig_width,
    #            margins[1] / fig_height,
    #            (fig_width - margins[2] - margins[0]) / fig_width,
    #            (fig_height - margins[3] - margins[1]) / fig_height
    #            ]

    ax.set_position(ax_rect)

    # Add axis titles and frame label; use absolute locations, rather then leaving it up to Matplotlib
    if ref is not None:
        plt.text(
            ax_rect[1] * 0.15,
            ax_rect[-1] + ax_rect[1],
            ref,
            horizontalalignment="left",
            verticalalignment="center",
            fontsize=fontsize * 1.2,
            transform=fig.transFigure,
        )
    plt.text(
        ax_rect[0] + 0.5 * ax_rect[2],
        0,
        xtitle,
        horizontalalignment="center",
        verticalalignment="bottom",
        fontsize=fontsize,
        transform=fig.transFigure,
    )
    plt.text(
        ax_rect[1] * 0.005,
        ax_rect[1] + 0.5 * ax_rect[3],
        ytitle,
        horizontalalignment="left",
        verticalalignment="center",
        fontsize=fontsize,
        transform=fig.transFigure,
        rotation=90,
    )

    if legend is not None:
        plt.rcParams["legend.frameon"] = False
        plt.rcParams["legend.edgecolor"] = "none"
        plt.rcParams["legend.labelspacing"] = 0.25
        plt.rcParams["legend.handlelength"] = 1
        plt.rcParams["legend.handletextpad"] = 0.25
        plt.legend(
            legend,
            prop={"size": fontsize * 0.8},
        )

    # fig.canvas.mpl_connect("key_press_event", key_press_action)

    if not boxed:
        ax.spines["right"].set_visible(False)
        ax.spines["top"].set_visible(False)

    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    ax.tick_params(
        which="both",
        axis="both",
        bottom=True,
        top=boxed,
        left=True,
        right=boxed,
        direction="in",
    )
    ax.tick_params(which="major", length=4, labelsize=fontsize * 0.8, width=1)
    ax.tick_params(which="minor", length=2, width=1)
    for axis in ["top", "bottom", "left", "right"]:
        ax.spines[axis].set_linewidth(1)

    if xlog:
        ax.semilogx()
    else:
        ax.ticklabel_format(axis="x", style="sci", scilimits=(-3, 3))
    if ylog:
        ax.semilogy()
    else:
        ax.ticklabel_format(axis="y", style="sci", scilimits=(-3, 3))

    try:
        ax.set_xlim([float(x) if x else None for x in xrange])
    except (ValueError, TypeError):
        pass
    try:
        ax.set_ylim([float(x) if x else None for x in yrange])
    except (ValueError, TypeError):
        pass

    if txt:
        pos = (ax.get_xlim()[1] * 1.025, ax.get_ylim()[1])
        if ylog:
            ypos = np.logspace(
                np.log10(ax.get_ylim()[1]), np.log10(ax.get_ylim()[0]), 25
            )
        else:
            ypos = np.linspace(ax.get_ylim()[1], ax.get_ylim()[0], 25)
        for t, y in zip(txt, ypos):
            plt.text(pos[0], y, t, size=fontsize / 1.5)

    if not GUI and save is None:
        plt.show()

    if save is not None:
        if not os.path.exists(os.path.dirname(save)):
            os.makedirs(os.path.dirname(save))
        base, ext = os.path.splitext(save)
        if ext == ".emf":
            save = base + ".pdf"
        fig.savefig(save, dpi=600, transparent=True)
        if ext == ".emf":
            try:
                subprocess.call(
                    [
                        "C:\Program Files\Inkscape\inkscape.exe",
                        "--file",
                        save,
                        "--export-emf",
                        base + ".emf",
                    ]
                )
                os.remove(save)
            except:
                print(
                    "Install Inkscape for conversion to emf.\nPlot saved as pdf in stead."
                )

    return fig, ax


def plot_heatmap(window, plot, xdata, ydata, size=100):
    ax = plt.gca()

    extent = ax.axis()
    plt.hexbin(xdata, ydata, cmap=plt.cm.Blues, vmax=3, vmin=0)
    return

    ax.xlim((-1, 3))
    if xrange is not None:
        xrange = ax.get_xlim()
    if yrange is not None:
        yrange = ax.get_ylim()

    print(xrange, yrange)


def change_extension(filename, new_extension, label=None, dir=None):
    if dir is not None:
        filename = f"{os.path.dirname(filename)}\\{dir}\\{os.path.basename(filename)}"
    if label is None:
        return f"{os.path.splitext(filename)[0]}.{new_extension}"
    else:
        return Path(f"{os.path.splitext(filename)[0]}_{label}.{new_extension}")


def create_new_path(base=r"c:/tmp/testPathlib/", folder="data", number=1):
    def new_path():
        return Path(
            base
            + os.getlogin()
            + "/"
            + datetime.date.today().strftime("%Y%m%d")
            + rf"/{folder}_{number:03}"
        )

    filepath = new_path()
    while filepath.exists():
        number += 1
        filepath = new_path()
    filepath.mkdir(parents=True, exist_ok=True)
    # print(filepath / 'test.txt')
    return filepath


def convert_to_h5(filename, check_existing=True, over_write=False):
    file_out = change_extension(filename, "h5")

    if check_existing or over_write:
        if os.path.isfile(file_out):
            if over_write:
                os.remove(file_out)
            else:
                print(f"{file_out} already exists, no conversion necessary")
                return file_out

    extension = os.path.splitext(filename)[1][1:]
    if extension == "dat":
        traces = pd.read_csv(filename, sep="\t", header=(0))
    elif extension == "tdms":
        tdms_file = TdmsFile(filename)
        tracked = tdms_file.as_dataframe()
        for i, trace in enumerate(tracked.items()):
            label = re.findall(r"\d+", trace[0])
            if label == []:
                label = " shared"
            else:
                label = label[0]
            channel = trace[0].replace(label, "")
            channel = channel.split("/")[-1][1:-1]
            h5_write_trace(trace[1].values, data, channel)
    else:
        print(f"Convert_to_h5: Format of file *.{extension} is not implemented!")

    # preliminary data processing
    _, labels, _ = h5_contents(file_out, to_file=True)
    for label in labels:
        h5_write_par(data, " selected", 1, type="local")

    if "Waad" in filename:
        z_multiplier = 1.333**2
        print(f"Correction factor {z_multiplier:.2f} applied to z-trace: {filename}")
    else:
        z_multiplier = 1.0

    for label in labels:
        for channel in ["X (um)", "Y (um)", "Z (um)"]:
            trace = h5_read_trace(data, channel)
            par = channel[0] + "0" + channel[1:]
            if channel == "Z (um)":
                h5_write_trace(trace * z_multiplier, data, channel)
                value = np.quantile(trace, 0.05)
                h5_write_par(
                    data, "dZ (um)", np.quantile(trace, 0.95) - value, type="local"
                )
            else:
                value = np.median(trace)
            h5_write_par(data, par, value, type="local")

    return file_out


def convert_log_h5(filename):
    """
    Convert init file to h5. Avoid characters '#' and '/', as these are not properly converted

    :param filename: File to be converted
    """
    f = configparser.ConfigParser()
    f.read(change_extension(filename, "log"))
    log = dict(f._sections)
    for k in log:
        log[k] = dict(f._defaults, **log[k])
        log[k].pop("__name__", None)

    hf = h5py.File(change_extension(filename, "h5"), "a")
    for section in log:
        for item in log[section]:
            entry = f"log/{section}/{item}"
            data = log[section][item]
            try:
                stored_data = hf[entry]
                stored_data[...] = data
            except KeyError:
                hf.create_dataset(entry, data=data)
    hf.close


# @pr.profile


def test(filename, label):
    # data = H5_data(filename, label)
    # data.save_xlsx()
    for label in tqdm(range(53)):
        data = H5_data(filename, label)
        data.save_xlsx()
    data.read_xlsx(36)
    print(data.parameters)
    print(data.traces)
    print(data.shared_traces)


if __name__ == "__main__":
    #   Intended as examples. Make your own script and import this module to your script.
    #   You can copy the code below as a start of your own script

    if False:
        filenames = [rf"C:\Users\noort\data\ForceTDMS\data_026.h5"]
        # filenames = [rf'C:\Users\noort\LTS\werk\Software\Python\testdata\ForceTDMS\\data_025.h5']
        for filename in filenames:
            data = Traces(filename=filename)
            # with data:            # for v in np.random.random(10):
            #     data.write_par('Z0 (um)', 0, 100 * v, 5 * v)
            # print(data.read_par(1, 'Z0 (um)'))
            # print(data.read_df('processed', filter='pN'))
            # t, l = data.contents()
            # print(t)
            # print(l)
            plt.plot(np.arange(0, 10, 1), np.random.random(10))
            # format_plot('x-as', 'y-ax', 'title', width=5, height=5)

    if False:
        filenames = [rf"C:\Users\noort\data\ForceDAT\data_002.dat"]

        for filename in filenames:
            filename = convert_dat_h5(filename)

            with Traces(filename) as data:
                for trace in data.contents():
                    print(trace)
                z = data.read_trace(
                    r"C:\Users\noort\data\ForceDAT\data_002.h5/data/processed/Z (um) > 91"
                )
                t = data.read_trace(
                    r"C:\Users\noort\data\ForceDAT\data_002.h5/data/processed/Time (s)"
                )
                plt.plot(t, z)
                format_plot("t (s)", "z (um)")

    # # file_contentsnew_files
    filename = r"data\data_010.tdms"
    filename = convert_to_h5(filename, over_write=True)
    # # channels, labels, parameters = h5_contents(filename, to_file=True)
    # # print(channels)
    # # print(labels)
    # # print(parameters)
    # x = h5_read_trace(filename, 'Time (s)')
    # print(x)
    filename = rf"data\data_014.h5"

    # test(filename, 18)
    # test(filename, 36)
    # test(filename, 19)
