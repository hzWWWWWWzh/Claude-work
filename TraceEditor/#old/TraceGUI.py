#!/usr/bin/env python
import copy
import inspect
import os
import re
import sys
import winsound
import warnings
from scipy.optimize import OptimizeWarning
from math import log10, floor
from pathlib import Path
from time import sleep

import PySimpleGUI as sg
import addcopyfighandler
import keyboard as kb
import numpy as np
import pandas as pd
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import matplotlib.pyplot as plt
import seaborn as sns

import ForceSpectroscopy as fs
import TraceIo as tio
import ctypes


# Tell Windows to scale the GUI properly on high resolution displays
ctypes.windll.user32.SetProcessDPIAware()

BG_COLOR = "white"
FG_COLOR = "lightgrey"
FG_COLOR = "#f6f6f6"

screen = {}
screen["size"] = sg.Window("").get_screen_size()
if screen["size"][0] > 1920:
    screen["plot_size"] = 0.9
    screen["plot_aspectratio"] = 0.83
    screen["scaling"] = 2.0
else:
    # print('no scaling')
    screen["plot_size"] = 0.95
    screen["plot_aspectratio"] = 0.75
    screen["scaling"] = 1.0


fitmodels = [
    "none",
    "ExponentialDrift",
    "MagnetForce",
    "zWLC",
    "zFiber",
    "zHMf",
    "zHMf2",
]


def handle_exceptions(func):
    def wrapper(*args, **kwargs):
        try:
            return func(*args, **kwargs)
        except Exception as e:
            exception_type, exception_object, exception_traceback = sys.exc_info()
            line_number = exception_traceback.tb_lineno
            print(">>> Exception type: ", exception_type)
            print(">>> File name: ", exception_traceback.tb_frame.f_code.co_filename)

    return wrapper


def generate_plot_text(window, data):
    if data.filename is not None:
        path = Path(str(window["file"].Get()))
        name = rf"{path.parts[-2]}/{path.parts[-1]}"
        name = name.split(".")[0].replace(" ", "_")

        text = rf'{name}: {window["label"].Get()}'
        length = 22  # split long filename
        text = [
            text[i * length : i * length + length]
            for i in range(1 + len(text) // length)
        ]

        par_list = data.list_parameters(format="fit_list")

        model_name = window["model"].get()
        if model_name != "none":
            _, _, names = tio.prepare_fit(model_name)
            for p in par_list:
                if p[0] == f" R2 {model_name}":
                    text.append(f"{model_name}: R2 = {p[1]}")

            text.append("")
            for p in par_list:
                if p[0] in list(names.values()):
                    text.append(f"{p[0]} = {p[1]}")
        else:
            text.append("")
            for p in par_list:
                if p[0][0] != " ":
                    text.append(f"{p[0]} = {p[1]}")
        return text
    else:
        return [""]


def get_index(items, item, default_index=0):
    if item in items:
        index = items.index(item)
    else:
        index = default_index
        index = np.clip(index, 0, len(items) - 1)
    return index


def get_default_pars(func):
    signature = inspect.signature(func)
    return {
        k: v.default
        for k, v in signature.parameters.items()
        if v.default is not inspect.Parameter.empty
    }


def adjust_value(value, increase=True, ctrl=False, shift=False):
    def extract_float(string):
        # The regular expression looks for any number, which may include a negative sign, may or may not include a decimal point
        match = re.search(r"[-+]?\d+(\.\d+)?", string)
        if match:
            return float(match.group(0))
        else:
            return None

    value = extract_float(value)
    if not value:
        return None

    if value == 0:
        # exponent = 0
        # value = 1
        return 1
    else:
        exponent = int(floor(log10(abs(value)))) - 1
    if increase:
        multiplier = 1
    else:
        multiplier = -1
    if ctrl:
        exponent -= 1
    if shift:
        exponent += 1
    value += multiplier * 10**exponent
    return value


def adjust_value_at_position(value_str, position, decrease=False):
    value_list = np.asarray(list(value_str))
    try:
        if value_list[position] == ".":
            position -= 1
    except IndexError:
        return value_str
    new_value_list = np.asarray(
        ["1" if i == position else "0" for i, v in enumerate(value_list)]
    )
    new_value_list[value_list == "."] = "."
    new_value_list[value_list == "-"] = "0"
    addition = "".join(new_value_list)
    if decrease:
        multiplier = -1
    else:
        multiplier = 1
    try:
        significance = len(value_str.split(".")[1])
    except IndexError:
        significance = 1
    return f"{float(value_str) + multiplier*float(addition):.{significance}f}"


def update_selection(window, plot, data, criteria):
    selected_all = None
    for channel, criterion in criteria.items():
        trace = data.traces[channel].values
        if channel == "Z (um)":
            trace -= data.calculate_drift()

        selected = np.ones_like(trace)
        if selected_all is None:
            selected_all = selected

        if criterion["der"]:
            trace = np.diff(trace)
            trace = np.append(trace, trace[-1])
        try:
            selected *= trace > float(criterion["min"])
        except ValueError:
            pass
        try:
            selected *= trace <= float(criterion["max"])
        except ValueError:
            pass

        selected_all *= selected
    if criterion["inv"]:
        selected_all = 1 - selected_all

    data.traces["Selection"] = selected_all
    update_plot(window, plot, data)
    return


@handle_exceptions
def computation_stack(
    window, plot=None, data=None, task=None, batch=False, parameters=None, stop=False
):
    # Use this thread for batch operations
    # print('@',task)
    # if (task is not None) and (parameters is not None):
    #     for par in parameters:
    #         print(par, parameters[par])

    # print('>>>>>>>>>>>>>>\n')
    # print(task)
    # print(batch)
    # print(parameters)
    # print('<<<<<<<<<<<<<<\n')

    global stack, stack_size, stack_params, results_df
    loop_time = 0.01

    if stop:
        stack = []

    ### fill stack with tasks ########################################################################################
    if task is not None:
        stack = []
        # print('>>>>>> Filling stack')
        window["progbar"].UpdateBar(100)
        window["computing"].Update("Computing " + task)
        results_df = None
        # fill stack
        if batch:
            stack_params = parameters
            for filename in window["file"].Widget["values"]:
                if task in ["MagnetForce2"]:
                    stack = stack.append([filename, "All", task])
                    # tio.debug([filename, 'All', task, stack])
                else:
                    labels = data.list_labels(selected_only=True)
                    if len(labels) == 0:
                        print(f"No selected labels found in file {filename}")
                    for label in labels:
                        stack.append([filename, label, task])
            stack_size = len(stack)
        else:
            stack_params = parameters
            stack.append([window["file"].Get(), window["label"].get(), task])
            stack_size = 1

        # print(f'Computing {stack_size} tasks: {task}')
        # for s in stack:
        #     print(s)
        return loop_time, data

    ### start of task processing ####################################################################################
    if len(stack) > 0:
        try:
            # process last item in stack
            # tio.debug(stack[-1], stack_params)
            filename, label, task = stack[-1]
            data.read(filename, label)

            if task == "MagnetForce":
                try:
                    shift = data.traces["Shift (mm)"]
                except KeyError:
                    print(f'No trace "Shift (mm)" found in file {filename}')
                    return loop_time
                model, _, names = tio.prepare_fit(task)
                for key in stack_params:
                    data.parameters.loc[names[key]] = stack_params[key].value
                data.traces["F (pN)"] = model.eval(stack_params, x=shift)

            elif task == "Select":
                if stack_params == ["All"]:
                    data.parameters[" selected"] = 1.0
                else:
                    data.parameters[" selected"] = 0.0
                data.parameters.loc["global", " selected"] = 0
                data.save(traces=False)
                update_GUI(window, data)

            elif task == "Select points":
                update_selection(window, plot, data, stack_params)
                data.save(parameters=False, traces=True)

            elif task == "Export":
                data.save_excel()

            elif task == "Select traces":
                selected = True
                for par, criteria in stack_params.items():
                    value = data.get_parameter(par)
                    select = criteria[0] < value <= criteria[1]
                    if criteria[-1]:
                        select = not select
                    selected = selected & select
                data.set_parameter(" selected", selected.astype(int), type="local")
                data.save(parameters=True, traces=False)

            elif task in fitmodels[1:]:
                ### fit model ################################################################################################
                model, fit_params, names = tio.prepare_fit(task, data)
                try:
                    x = data.traces[names["x"]].values.copy()
                    y = data.traces[names["y"]].values.copy()
                except KeyError:
                    print(
                        f'No trace "{names["x"]}" or "{names["y"]}" found in file {filename}'
                    )
                    return loop_time

                if task == "ExponentialDrift":
                    # Only use data at low force and large enough amplitude
                    amplitude = data.traces["A (a.u.)"]
                    time = data.traces["Time (s)"]
                    ruptured = np.asarray(
                        amplitude <= np.percentile(amplitude, 90) / 10
                    )
                    if ruptured.any():
                        selection = time < (np.min(time[ruptured])) - 2.0
                    else:
                        selection = time > 0
                    if "Shift (mm)" in data.list_channels():
                        shift = data.traces["Shift (mm)"]
                        selection = selection & (shift >= np.percentile(shift, 90))
                    if "Rotation (turns)" in data.list_channels():
                        rotation = np.abs(data.traces["Rotation (turns)"])
                        selection = selection & (
                            rotation >= np.percentile(rotation, 90)
                        )

                    data.traces["Selection"] = selection
                    data.set_parameter(
                        "Z0 (um)",
                        np.median(data.traces["Z (um)"]),
                        type="local",
                        error=None,
                    )
                    stack_params["Z0"]["value"] = data.get_parameter("Z0 (um)")
                else:
                    if names["y"] == "Z (um)":
                        y -= data.calculate_drift()
                        if task != "zWLC":
                            y += data.get_parameter("Zbead (um)")
                    if names["x"] == "Z (um)":
                        x -= data.calculate_drift() + data.get_parameter("Zbead (um)")
                    selection = data.traces["Selection"].values.astype(bool)

                for key, p in stack_params.items():
                    if p["fit"] == "fit":
                        fit_params[key].value = p["value"]
                        fit_params[key].vary = True
                        fit_params[key].min = p["min"]
                        fit_params[key].max = p["max"]
                    elif p["fit"] == "local":
                        fit_params[key].value = data.get_parameter(
                            names[key], default=p["value"]
                        )
                        fit_params[key].vary = False
                    else:
                        fit_params[key].value = p["value"]
                        fit_params[key].vary = False

                if (
                    np.any([fit_params[key].vary for key in fit_params])
                    and selection.any()
                ):
                    y = np.asarray(y)
                    with warnings.catch_warnings():
                        warnings.simplefilter("ignore", category=OptimizeWarning)
                        result = model.fit(y[selection], fit_params, x=x[selection])

                    for key, p in stack_params.items():
                        error = result.params[key].stderr if p["fit"] == "fit" else None
                        data.set_parameter(
                            names[key],
                            result.params[key].value,
                            error=error,
                            type=p["fit"],
                        )

                    if task == "ExponentialDrift":
                        if ruptured.any():
                            selection = time < (np.min(time[ruptured])) - 2.0
                        else:
                            selection = time > 0
                        trace = data.traces["Z (um)"][selection] - model.eval(
                            result.params, x=data.traces["Time (s)"][selection]
                        )
                        offset = np.percentile(trace, 2)
                        data.set_parameter(
                            "Z0 (um)", result.params["Z0"] + offset, type="local"
                        )

                        # Calculate offsets in x and y
                        for axe in ["X", "Y"]:
                            mean = np.median(data.traces[f"{axe} (um)"][selection])
                            std = np.std(data.traces[f"{axe} (um)"][selection])
                            data.set_parameter(
                                f"{axe}0 (um)", mean, type="fit", error=std
                            )

                        # Calculate dZ
                        selection = amplitude >= np.percentile(amplitude, 90) / 100
                        trace = data.traces["Z (um)"] - data.calculate_drift()
                        dz = np.percentile(trace, 98)
                        data.set_parameter("dZ (um)", dz, type="local")

                        try:
                            r2 = 1 - result.residual.var() / np.var(y[selection])
                            n = len(x[selection])
                            p = len(result.covar) - 1
                            r2_adj = r2 - (1 - r2) * p / (n - p - 1)
                            data.set_parameter(f" R2 {task}", r2_adj)
                        except TypeError:
                            pass
                else:
                    for key, p in stack_params.items():
                        data.set_parameter(
                            names[key], fit_params[key].value, type=p["fit"]
                        )
                data.save(traces=False)

            if task == "Delete":
                pars = [
                    c
                    for c in data.list_parameters(format="names")
                    if c not in stack_params
                ]
                data.parameters = data.parameters[pars]
                data.save(parameters=True, traces=False)

            if task == "Phase jump":
                x = tio.h5_read_trace(data, "X (um)")
                y = tio.h5_read_trace(data, "Y (um)")
                z = tio.h5_read_trace(data, "Z (um)")
                tio.h5_write_trace(fs.correct_phase_jump(x, y, z), data, "Z (um)")

            if task == "Free DNA":
                f = data.traces["F (pN)"].values
                z = (
                    data.traces["Z (um)"].values
                    - data.calculate_drift()
                    - data.get_parameter("Zbead (um)")
                )
                L = stack_params["L (bp)"]
                A = stack_params["A (nm)"]
                S = stack_params["S (pN)"]

                data.traces["Lfree (bp)"] = fs.get_free_DNA(z, f, L, A, S)
                data.save(parameters=False, traces=True)

            if task == "Z bead":
                if "Shift (mm)" in data.list_channels():
                    f = data.traces["F (pN)"].values
                    z = data.traces["Z (um)"].values - data.calculate_drift()
                    L = stack_params["L (bp)"]
                    A = stack_params["A (nm)"]
                    S = stack_params["S (pN)"]

                    amplitude = data.traces["A (a.u.)"].values
                    time = data.traces["Time (s)"]
                    ruptured = np.asarray(
                        amplitude <= np.percentile(amplitude, 90) / 10
                    )
                    if ruptured.any():
                        selection = time < (np.min(time[ruptured])) - 2.0
                    else:
                        selection = time > 0

                    fit = fs.fit_zbead(f, z, L, A, S, selection=selection)
                    if fit[1] < 0.5:
                        data.set_parameter("Zbead (um)", *fit, type="fit")
                        for name, value in stack_params.items():
                            data.set_parameter(name, value, type="global")
                        data.save(parameters=True, traces=False)
                if "Rotation (turns)" in data.list_channels():
                    data = fs.fit_zbead_from_rotation(data)
                    data.save(parameters=True, traces=False)

            if task == "Ruptures":
                df = fs.get_ruptures(filename, label)
                try:
                    results_df = results_df.append(df)
                except AttributeError:
                    results_df = df
            if task == "Heatmap":
                cs = plot["ax"].collections[1]
                data = np.asarray(cs.get_offsets())
                df = pd.DataFrame(data, columns=["xdata", "ydata"])
                try:
                    results_df = results_df.append(df, ignore_index=True)
                except AttributeError:
                    results_df = df

            if task == "Save graph":
                if stack_params is None:
                    fileout = filename.split(".")[0] + "_{label}.jpg"
                    fileout = sg.popup_get_file(
                        "Give filename template",
                        save_as=True,
                        no_window=True,
                        default_path=fileout,
                        background_color=FG_COLOR,
                    )
                    fileout = os.path.abspath(fileout).replace(
                        filename.split(".")[0], "{filename}"
                    )
                    stack_params = fileout

                fileout = stack_params.replace("{label}", label)
                fileout = fileout.replace("{filename}", filename.split(".")[0])

                if Path(fileout).is_file():
                    fileout = sg.popup_get_file(
                        "Overwrite file?",
                        save_as=True,
                        no_window=True,
                        default_path=fileout,
                        background_color=FG_COLOR,
                    )
        except Exception as e:
            exception_type, exception_object, exception_traceback = sys.exc_info()
            line_number = exception_traceback.tb_lineno
            print(">>> Exception type: ", exception_type)
            print(">>> File name: ", exception_traceback.tb_frame.f_code.co_filename)
            print(">>> Line number: ", line_number)

        ### update GUI items ##########################################################################################

        window["progbar"].UpdateBar(100 * len(stack) / stack_size)
        update_GUI(window, data)
        update_plot(window, plot, data)
        stack = np.delete(stack, -1, axis=0)
        return loop_time, data
    else:
        window["computing"].Update("Idle")
        window["progbar"].UpdateBar(0)
        winsound.MessageBeep()
        return None, data


@handle_exceptions
def GUI_select_points(window, data, plot):
    filename = window["file"].Get()
    label = window["label"].get()
    data.read(filename, label)
    channels = data.list_channels()

    if "Selection" in channels:
        old_selection = data.traces["Selection"].values.astype(bool).copy()
    else:
        old_selection = np.ones_like(data.traces["Z (um)"].values)

    def_size = (10, 0)
    layout = [
        [
            sg.Text(
                " " * 27 + "min" + " " * 24 + "max" + " " * 12 + "derivative  invert",
                pad=(10, 5),
            )
        ]
    ]
    for channel in channels:
        layout.append(
            [
                sg.Text(channel, pad=(10, 5), size=def_size),
                sg.InputText(
                    "",
                    key="min_" + channel,
                    pad=(15, 5),
                    background_color=BG_COLOR,
                    enable_events=True,
                    size=(12, 0),
                ),
                sg.InputText(
                    "",
                    key="max_" + channel,
                    pad=(15, 5),
                    background_color=BG_COLOR,
                    enable_events=True,
                    size=(12, 0),
                ),
                sg.Checkbox(
                    "",
                    default=False,
                    key="der_" + channel,
                    size=(0, 10),
                    enable_events=True,
                    background_color=FG_COLOR,
                    text_color="black",
                ),
                sg.Checkbox(
                    "",
                    default=False,
                    key="inv_" + channel,
                    size=(0, 10),
                    enable_events=True,
                    background_color=FG_COLOR,
                    text_color="black",
                ),
            ]
        )

    layout.append(
        [
            sg.Button(
                "Cancel", button_color=("black", FG_COLOR), size=def_size, pad=(10, 5)
            ),
            sg.Button(
                "Batch", button_color=("black", FG_COLOR), size=def_size, pad=(10, 5)
            ),
            sg.Button("OK", pad=(10, 5), size=def_size),
        ]
    )

    selection_window = sg.Window(
        "Select range",
        layout,
        default_element_size=(10, 0),
        border_depth=1,
        return_keyboard_events=True,
        location=(925, 20),
        finalize=True,
    )

    selection_data = {}
    criteria = {}

    for channel in channels:
        selection_data[channel] = data.traces[channel].values.copy()
        if channel == "Z (um)":
            selection_data[channel] -= data.calculate_drift() + data.get_parameter(
                "Zbead (um)"
            )

    while True:
        event, values = selection_window.Read()
        for channel in channels:
            criteria[channel] = {
                "min": values["min_" + channel],
                "max": values["max_" + channel],
                "der": values["der_" + channel],
                "inv": values["inv_" + channel],
            }
        criteria = {
            k: v
            for k, v in criteria.items()
            if v != {"min": "", "max": "", "der": False, "inv": False}
        }

        try:
            focus = selection_window.FindElementWithFocus().Key
        except:
            focus = None

        # if event in ['Left:37', 'Right:39']:
        #     value = selection_window[focus].get()
        #     try:
        #         value = adjust_value(value, 'Right' in event, kb.is_pressed('ctrl'), kb.is_pressed('shift'))
        #         value = f'{value:.4g}'
        #         selection_window[focus].update(value=value)
        #     except:
        #         pass

        if event[:4] in ["min_", "max_", "der_", "inv_"]:
            update_selection(window, plot, data, criteria)

        if event in ["Batch", "OK"]:
            selection_window.close()
            return criteria, event == "Batch"

        elif event == "Cancel":
            # tio.h5_write_trace(old_selection, data, 'Selection')
            data.traces["Selection"] = old_selection
            selection_window.close()
            update_plot(window, plot, data)
            return None, None


def GUI_fit_parameters(window, data, plot):
    if window["model"].get() == "none":
        print("No model selected")
        return None, None

    old_parameters = data.parameters.copy()
    previous_parameters = data.parameters.copy()

    model_name = window["model"].get()
    _, model_params, names = tio.prepare_fit(model_name, data)

    sg.theme("LightGrey2")
    def_size = (10, 0)
    layout = []
    type_options = ["global", "local", "fit"]
    layout.append(
        [
            sg.Text("", size=def_size),
            sg.Text("Set all to  ", size=(9, 0)),
            sg.Combo(
                type_options,
                default_value="choose",
                background_color=BG_COLOR,
                key="reset",
                enable_events=True,
                size=def_size,
            ),
            sg.Text("min", size=(9, 0)),
            sg.Text("max", size=(9, 0)),
        ]
    )
    for i, key in enumerate(model_params):
        try:
            type = (
                "fit"
                if pd.notna(data.parameters.loc[data.label, "e:" + names[key]])
                else "local"
            )
        except KeyError:
            type = "local"

        layout.append(
            [
                sg.Text(names[key], size=def_size),
                sg.InputText(
                    model_params[key].value,
                    key="val_" + key,
                    background_color=BG_COLOR,
                    enable_events=True,
                    size=def_size,
                ),
                sg.Combo(
                    type_options,
                    default_value=type,
                    background_color=BG_COLOR,
                    key="fit_" + key,
                    enable_events=True,
                    size=def_size,
                ),
                sg.InputText(
                    model_params[key].min,
                    key="min_" + key,
                    background_color=BG_COLOR,
                    enable_events=True,
                    size=def_size,
                ),
                sg.InputText(
                    model_params[key].max,
                    key="max_" + key,
                    background_color=BG_COLOR,
                    enable_events=True,
                    size=def_size,
                ),
            ]
        )

    layout.append(
        [
            sg.Button("Cancel", button_color=("black", "lightgrey"), size=def_size),
            sg.Button("Batch", button_color=("black", "lightgrey"), size=def_size),
            sg.Button("OK", button_color=("black", "lightgrey"), size=def_size),
        ]
    )

    fit_window = sg.Window(
        "Fit parameters " + model_name,
        layout,
        default_element_size=(10, 0),
        border_depth=1,
        return_keyboard_events=True,
        location=(925, 20),
        finalize=True,
    )

    while True:
        event, values = fit_window.Read()
        try:
            focus = fit_window.FindElementWithFocus().Key
        except:
            focus = None
        if event is None:
            break

        elif event == "reset":
            for key in model_params:
                fit_window[f"fit_{key}"].update(
                    set_to_index=type_options.index(fit_window["reset"].get())
                )

        elif (
            "val_" in focus
            and event not in ["Cancel", "Batch", "OK"]
            and focus != event
        ):
            if event in ["Up:38", "Down:40"]:
                cursor_pos = (
                    fit_window[focus].Widget.index("insert")
                    if fit_window[focus].Widget
                    else -1
                )
                value = adjust_value_at_position(
                    fit_window[focus].get(), cursor_pos - 1, "Down" in event
                )
                fit_window[focus].update(value=value)
                fit_window[focus].Widget.icursor(cursor_pos)

        elif event in ["OK", "Batch"]:
            fit_params = {}
            for key in sorted(model_params):
                fit_params[key] = {
                    "value": float(fit_window["val_" + key].get()),
                    "fit": fit_window["fit_" + key].get(),
                    "min": float(fit_window["min_" + key].get()),
                    "max": float(fit_window["max_" + key].get()),
                }

            fit_window.close()
            return fit_params, event == "Batch"

        elif event == "Cancel":
            data.parameters = old_parameters
            fit_window.close()
            update_plot(window, plot, data)
            return None, None

        for key in model_params:
            try:
                data.set_parameter(names[key], float(fit_window["val_" + key].get()))
            except ValueError:
                pass

        if not previous_parameters.equals(data.parameters):
            update_plot(window, plot, data)
            previous_parameters = data.parameters.copy()


def GUI_fix_parameters(window, data, plot, params=None):

    if params is None:
        return None, None

    old_parameters = data.parameters.copy()

    sg.theme("LightGrey2")
    def_size = (10, 0)
    layout = []

    for name, value in params.items():
        value = data.get_parameter(name, default=value)
        layout.append(
            [
                sg.Text(name, size=def_size),
                sg.InputText(
                    value,
                    key="val_" + name,
                    background_color=BG_COLOR,
                    enable_events=True,
                    size=def_size,
                ),
            ]
        )

    layout.append(
        [
            sg.Button("Cancel", button_color=("black", "lightgrey"), size=def_size),
            sg.Button("Batch", button_color=("black", "lightgrey"), size=def_size),
            sg.Button("OK", button_color=("black", "lightgrey"), size=def_size),
        ]
    )

    fit_window = sg.Window(
        "Fix parameters ",
        layout,
        default_element_size=(10, 0),
        border_depth=1,
        return_keyboard_events=True,
        location=(925, 20),
        finalize=True,
    )

    while True:
        event, values = fit_window.Read()
        try:
            focus = fit_window.FindElementWithFocus().Key
        except:
            focus = None
        if event is None:
            break

        elif event in ["OK", "Batch"]:
            for name, value in params.items():
                params[name] = float(fit_window["val_" + name].get())
            fit_window.close()
            return params, event == "Batch"

        elif event == "Cancel":
            data.parameters = old_parameters
            fit_window.close()
            return None, None


def GUI_select_traces(window, data):
    parameters = data.list_parameters(format="names")
    if " selected" in parameters:
        parameters.remove(" selected")

    def_size = (10, 0)
    layout = [
        [
            sg.Text(
                " " * 27 + "min" + " " * 24 + "max" + " " * 12 + "      invert",
                pad=(10, 5),
            )
        ]
    ]
    for par in parameters:
        layout.append(
            [
                sg.Text(par, pad=(10, 5), size=def_size),
                sg.InputText(
                    "",
                    key="min_" + par,
                    pad=(15, 5),
                    background_color=BG_COLOR,
                    enable_events=True,
                    size=(12, 0),
                ),
                sg.InputText(
                    "",
                    key="max_" + par,
                    pad=(15, 5),
                    background_color=BG_COLOR,
                    enable_events=True,
                    size=(12, 0),
                ),
                sg.Checkbox(
                    "",
                    default=False,
                    key="inv_" + par,
                    size=(0, 10),
                    enable_events=True,
                    background_color=FG_COLOR,
                    text_color="black",
                ),
            ]
        )

    layout.append(
        [
            sg.Button("Cancel", button_color=FG_COLOR, size=def_size, pad=(10, 5)),
            sg.Text(" ", size=(12, 0), pad=(10, 5)),
            sg.Button("OK", pad=(10, 5), size=def_size),
        ]
    )

    selection_window = sg.Window(
        "Select criteria",
        layout,
        default_element_size=(10, 0),
        border_depth=1,
        return_keyboard_events=True,
        location=(925, 20),
        finalize=True,
    )

    while True:
        event, values = selection_window.Read()
        try:
            focus = selection_window.FindElementWithFocus().Key
        except:
            focus = None

        if event in ["Left:37", "Right:39"]:
            value = selection_window[focus].get()

            value = adjust_value(
                value, "Right" in event, kb.is_pressed("ctrl"), kb.is_pressed("shift")
            )
            if value:
                value = f"{value:.4g}"
                selection_window[focus].update(value=value)

        if event in ["OK"]:
            criteria = {}
            default = [-np.inf, np.inf, 0]
            for par in parameters:
                selection_range = copy.copy(default)
                for i, c in enumerate(["min_", "max_"]):
                    try:
                        selection_range[i] = float(values[c + par])
                    except ValueError:
                        pass
                if values["inv_" + par]:
                    selection_range[-1] = 1
                criteria[par] = selection_range
            selection_window.close()
            criteria = {k: v for k, v in criteria.items() if v != default}
            return criteria

        elif event == "Cancel":
            selection_window.close()
            return


def GUI_histogram(window, plot=None, x=None, selected=None):

    if x is None:
        return

    sg.theme("LightGrey2")
    def_size = (10, 0)

    layout = []
    entrees = ["min", "max", "#bins", "N max"]

    try:
        hist_all, bins = np.histogram(x)
    except ValueError:
        return
    values = [bins[0], bins[-1], len(hist_all), np.max(hist_all) * 1.2]

    for entree, value in zip(entrees, values):
        layout.append(
            [
                sg.Text(entree, size=def_size),
                sg.InputText(
                    f"{value:.1f}",
                    key=entree,
                    background_color=BG_COLOR,
                    enable_events=True,
                    size=(12, 0),
                ),
            ]
        )
    layout.append(
        sg.Button("Update", button_color=("black", "lightgrey"), size=def_size)
    )

    hist_window = sg.Window(
        "Histogram",
        [layout],
        default_element_size=(20, 0),
        border_depth=1,
        return_keyboard_events=True,
        location=(925, 20),
        finalize=True,
    )
    while True:
        event, values = hist_window.Read()
        if event == "Update":
            bins = np.linspace(
                float(hist_window["min"].get()),
                float(hist_window["max"].get()),
                int(float(hist_window["#bins"].get())) + 1,
            )

            hist_all, bins = np.histogram(x, bins=bins)
            hist, bins = np.histogram(x[selected.astype(bool)], bins=bins)

            files = window["file"].Widget["values"]
            text = ["Files:"]
            for file in files:
                path = Path(file)
                name = rf"{path.parts[-2]}/{path.parts[-1]}"
                name = name.split(".")[0].replace(" ", "_")
                text.append(name)
            text.append("")
            text.append(f"N = {len(x[selected.astype(bool)])}")
            text.append(f"mean = {np.mean(x[selected.astype(bool)]):.2f}")
            text.append(f"std = {np.std(x[selected.astype(bool)]):.2f}")

            plot["ax"].cla()
            plot["ax"].stairs(hist_all, bins, fill=True, color="lightgray")
            plot["ax"].stairs(hist, bins, fill=True, color=window["color"].get())

            tio.format_plot(
                x.name,
                "N (-)",
                fig=plot["fig"],
                txt=text,
                ax=plot["ax"],
                xrange=[bins[0], bins[-1]],
                yrange=[0, float(hist_window["N max"].get())],
                GUI=True,
                scale_page=0.9,
                aspect=0.83,
            )
            plot["canvas"].draw()
            # plt.suptitle(f'N = {len(df[df[" selected"] == 0])}/{len(df)}', ha='left', x=0.02)

        elif event == sg.WIN_CLOSED:
            break


def update_plot(window, plot, data, save=False):
    if not Path(window["file"].get()).exists() or data.open_files == []:
        print(f'File {window["file"].get()} does not exist')
        plot["ax"].clear()
        plot["canvas"].draw()
        return

    if window["channel_x"].get() != "" and window["channel_y"].get() != "":
        plot_names, plot_data, plot_ranges = [], [], []
        for axe in ["x", "y"]:
            plot_names.append(window[f"channel_{axe}"].get())
            plot_data.append(data.traces[plot_names[-1]].copy())
            if (
                plot_names[-1] == "Z (um)"
                and window["model"].get() != "ExponentialDrift"
            ):
                plot_data[-1] -= data.calculate_drift() + data.get_parameter(
                    "Zbead (um)", default=0
                )
            elif plot_names[-1] != "Z (um)":
                offset = plot_names[-1][0] + "0" + plot_names[-1][1:]
                if offset in data.list_parameters(format="names"):
                    plot_data[-1] -= data.get_parameter(offset)
            plot_ranges.append(
                [window[f"channel_{axe}_min"].Get(), window[f"channel_{axe}_max"].Get()]
            )
        try:
            selected = data.traces["Selection"] == 1
        except KeyError:
            # selected = np.ones_like(data.traces.index) == 1
            selected = np.zeros_like(data.traces.index)

        # plot data
        plot["ax"].cla()
        plot["ax"].scatter(*plot_data, s=20, facecolors="none", edgecolors="lightgrey")
        plot["ax"].scatter(
            plot_data[0][selected],
            plot_data[1][selected],
            s=20,
            facecolors="none",
            edgecolors=window["color"].get(),
        )

        # plot model
        model_name = window["model"].get()

        if model_name not in ["none"]:
            model, params, names = tio.prepare_fit(model_name)
            if names["y"] in plot_names and names["x"] in data.traces.columns:
                for p in params:
                    if names[p] in data.list_parameters(format="names"):
                        params[p].value = data.get_parameter(names[p])
                if model_name == "zWLC":
                    params["z_bead"].value = 0
                model_data = model.eval(x=data.traces[names["x"]], params=params)

                model_plot = plot_data.copy()
                model_plot[plot_names.index(names["y"])] = model_data
                plot["ax"].plot(*model_plot, color="black")

        # format plot
        tio.format_plot(
            *plot_names,
            fig=plot["fig"],
            ax=plot["ax"],
            xlog=window["xlog"].get(),
            ylog=window["ylog"].get(),
            txt=generate_plot_text(window, data),
            xrange=plot_ranges[0],
            yrange=plot_ranges[1],
            GUI=True,
            scale_page=screen["plot_size"],
            aspect=screen["plot_aspectratio"],
        )
        plot["canvas"].draw()


def update_GUI(window, data, pars_only=False):
    if len(data.open_files) == 0:
        for element in ["file", "label", "params", "channel_x", "channel_y"]:
            window[element].update(values=[])
        return
    window["file"].update(
        values=data.open_files, set_to_index=data.open_files.index(data.filename)
    )

    if not pars_only:
        labels = data.list_labels(selected_only=window[" selected"].get())
        index = labels.index(data.label) if data.label in labels else 0
        window["label"].update(values=labels, set_to_index=index)

    # models
    model_name = window["model"].get()
    if model_name == "none":
        selection = [p for p in data.list_parameters(format="names") if " R2" not in p]
    else:
        _, _, names = tio.prepare_fit(model_name)

        parameter_list = data.list_parameters(format="names")
        selection = [" selected", " R2 " + model_name]
        selection += [item for item in names.values() if item in parameter_list]

    window["params"].update(
        values=data.list_parameters(format="list", selection=selection)
    )

    # channels
    channels = data.list_channels()
    for channel in ["channel_x", "channel_y"]:
        try:
            index = channels.index(window[channel].get())
        except ValueError:
            index = 0
        window[channel].update(values=channels, set_to_index=index)


def trace_editor(start_file=None, start_label=None):

    data = tio.hdf_data()
    col_width = 15
    pad = (5, 2)

    menu_def = [
        [
            "&File",
            [
                "&Open",
                "Import",
                "Close",
                ["Current file", "All files"],
                "&Save",
                ["Save all", "Export"],
                "Display h5",
                "E&xit",
            ],
        ],
        ["Select", ["Traces", ["All", "None", "Toggle", "Criteria..."], "Points..."]],
        ["Parameters", ["Adjust...", "Plot", "Delete"]],
        ["Process", ["Z bead", "Free DNA", "Phase jump", "Ruptures"]],
        [
            "&Plot",
            [
                "&Copy",
                "Type",
                ["Z(t)", "F(Z)", "Y(X)"],
                "Scale",
                ["Auto", "Fix", "Contour length"],
                "Play",
                "Heatmap",
            ],
        ],
        ["&Help", ["Break", "&About...", "Reduce Scaling", "Test"]],
    ]

    parameter_table = sg.Table(
        values=[],
        enable_events=True,
        headings=["parameter", "value"],
        display_row_numbers=False,
        auto_size_columns=False,
        col_widths=[11, 22],
        num_rows=25,
        size=(40, 100),
        key="params",
        justification="left",
        background_color=BG_COLOR,
        text_color="black",
        header_background_color=FG_COLOR,
        bind_return_key=True,
        change_submits=False,
    )

    parameter_column = [
        [
            sg.Combo(
                list(fitmodels),
                key="model",
                default_value=list(fitmodels)[0],
                size=[40, 10],
                pad=pad,
                enable_events=True,
                background_color=BG_COLOR,
                text_color="black",
            )
        ],
        [sg.Text("", pad=pad, background_color=BG_COLOR)],
        [parameter_table],
    ]

    filenames = sg.Combo(
        values=[],
        enable_events=True,
        background_color=BG_COLOR,
        size=[85, 35],
        pad=pad,
        key="file",
    )

    graph = sg.Canvas(size=(600, 500), key="canvas", background_color=BG_COLOR)
    graph_column = [
        [filenames],
        [graph],
        [
            sg.ProgressBar(
                100,
                size=(48 * screen["scaling"], 5 * screen["scaling"]),
                key="progbar",
                bar_color=("blue", FG_COLOR),
            )
        ],
        [
            sg.InputText(
                "Starting up",
                key="computing",
                size=[80, 5],
                background_color=BG_COLOR,
                text_color="black",
                border_width=0,
                visible=True,
            )
        ],
    ]

    plot_column = [
        [
            sg.Combo(
                [],
                key="label",
                background_color=BG_COLOR,
                size=[col_width, 10],
                pad=pad,
                enable_events=True,
            )
        ],
        [
            sg.Checkbox(
                " selected",
                default=False,
                key=" selected",
                size=(0, 10),
                pad=pad,
                enable_events=True,
                background_color=BG_COLOR,
                text_color="black",
            )
        ],
        [
            sg.Text(
                "___________", pad=pad, background_color=BG_COLOR, text_color="black"
            )
        ],
        [sg.Text("x-axis", pad=pad, background_color=BG_COLOR, text_color="black")],
        [
            sg.Combo(
                [""],
                key="channel_x",
                background_color=BG_COLOR,
                size=[col_width, 10],
                pad=pad,
                enable_events=True,
            )
        ],
        [sg.Text("min", pad=pad, background_color=BG_COLOR, text_color="black")],
        [
            sg.InputText(
                "",
                key="channel_x_min",
                size=[col_width, 10],
                pad=pad,
                enable_events=True,
                background_color=BG_COLOR,
                text_color="black",
            )
        ],
        [sg.Text("max", pad=pad, background_color=BG_COLOR, text_color="black")],
        [
            sg.InputText(
                "",
                key="channel_x_max",
                size=[col_width, 10],
                pad=pad,
                enable_events=True,
                background_color=BG_COLOR,
                text_color="black",
            )
        ],
        [
            sg.Checkbox(
                "Log",
                default=False,
                key="xlog",
                size=(0, 10),
                pad=pad,
                enable_events=True,
                background_color=BG_COLOR,
                text_color="black",
            )
        ],
        [
            sg.Text(
                "___________", pad=pad, background_color=BG_COLOR, text_color="black"
            )
        ],
        [sg.Text("y-axis", pad=pad, background_color=BG_COLOR, text_color="black")],
        [
            sg.Combo(
                [""],
                key="channel_y",
                size=[col_width, 10],
                pad=pad,
                enable_events=True,
                background_color=BG_COLOR,
                text_color="black",
            )
        ],
        [sg.Text("min", pad=pad, background_color=BG_COLOR, text_color="black")],
        [
            sg.InputText(
                "",
                key="channel_y_min",
                size=[col_width, 10],
                pad=pad,
                enable_events=True,
                background_color=BG_COLOR,
                text_color="black",
            )
        ],
        [sg.Text("max", pad=pad, background_color=BG_COLOR, text_color="black")],
        [
            sg.InputText(
                "",
                key="channel_y_max",
                size=[col_width, 10],
                pad=pad,
                enable_events=True,
                background_color=BG_COLOR,
                text_color="black",
            )
        ],
        [
            sg.Checkbox(
                "Log",
                default=False,
                key="ylog",
                size=(0, 10),
                pad=pad,
                enable_events=True,
                background_color=BG_COLOR,
                text_color="black",
            )
        ],
        [
            sg.Text(
                "___________________",
                pad=pad,
                background_color=BG_COLOR,
                text_color="black",
            )
        ],
        [sg.Text("Color", pad=pad, background_color=BG_COLOR, text_color="black")],
        [
            sg.Combo(
                ["red", "blue", "green", "orange", "purple"],
                default_value="blue",
                key="color",
                size=[col_width, 10],
                pad=pad,
                enable_events=True,
                background_color=BG_COLOR,
                text_color="black",
            )
        ],
        # [sg.Checkbox('Selection mode', default=True, key='selection mode', size=(0, 10), pad=pad, enable_events=True,
        #              background_color=BG_COLOR, text_color='black')],
        [],
    ]

    column_height = 570 * screen["scaling"]
    layout = [
        [sg.Menu(menu_def, tearoff=False)],
        [
            sg.Column(
                parameter_column,
                pad=pad,
                background_color=BG_COLOR,
                size=(310 * screen["scaling"], column_height),
            ),
            sg.Column(
                graph_column,
                pad=pad,
                background_color=BG_COLOR,
                size=(650 * screen["scaling"], column_height),
            ),
            sg.Column(
                plot_column,
                pad=pad,
                background_color=BG_COLOR,
                size=(130 * screen["scaling"], column_height),
            ),
        ],
    ]

    sg.theme("LightGrey2")
    window = sg.Window(
        "TraceEditor",
        layout,
        default_element_size=(20, 5),
        element_padding=pad,
        border_depth=1,
        return_keyboard_events=True,
        location=(0, 10),
        finalize=True,
        background_color=FG_COLOR,
    )

    plot = {}
    plot["window"] = window["canvas"]
    plot["fig"], plot["ax"] = tio.format_plot(empty_canvas=True)
    plot["canvas"] = FigureCanvasTkAgg(plot["fig"], plot["window"].TKCanvas)
    plot["canvas"].get_tk_widget().pack(side="top", fill="both", expand=1)
    plot["canvas"].draw()

    def on_click(event):
        global start_pos
        start_pos = (event.xdata, event.ydata)

    def on_release(event):
        if event.button == 3:
            update_plot(window, plot, data)
        else:
            end_pos = (event.xdata, event.ydata)
            if start_pos != end_pos:
                try:
                    plot["ax"].set_xlim(sorted([start_pos[0], end_pos[0]]))
                    plot["ax"].set_ylim(sorted([start_pos[1], end_pos[1]]))
                    plot["fig"].canvas.draw_idle()
                except TypeError:
                    pass

    plot["fig"].canvas.mpl_connect("button_press_event", on_click)
    plot["fig"].canvas.mpl_connect("button_release_event", on_release)

    update_GUI(window, data)
    update_plot(window, plot, data)
    window["computing"].Update("Idle")

    if start_file is not None and os.path.isfile(start_file):
        if data.read(start_file, start_label):
            window["channel_x"].update(
                values=data.list_channels(),
                set_to_index=data.list_channels().index("Time (s)"),
            )
            window["channel_y"].update(
                values=data.list_channels(),
                set_to_index=data.list_channels().index("Z (um)"),
            )
            update_GUI(window, data)
            update_plot(window, plot, data)

    wait_time = None
    while True:
        # try:
        event, values = window.Read(timeout=wait_time, timeout_key="Timeout")
        try:
            focus = window.FindElementWithFocus().Key
        except:
            focus = None

        # print(f'Event: {event}, focus: {focus}')

        if event == sg.WIN_CLOSED or event == "Exit" or event == "__TIMEOUT__":
            break

        elif event == "Timeout":
            wait_time, data = computation_stack(window, plot, data)

        # Menu events
        elif event in ["Open", "Import"]:
            extensions = ["*.hdf"] if event == "Open" else ["*.dat", "*.tdms", "*.h5"]

            filenames = list(
                sg.popup_get_file(
                    "filename to open",
                    no_window=True,
                    multiple_files=True,
                    file_types=(("trace files", extensions),),
                )
            )
            if len(filenames) > 0:
                window["progbar"].UpdateBar(100)
                window["computing"].Update("Opening files ...")
                for i, filename in enumerate(filenames):
                    file_type = Path(filename).suffix
                    if file_type == ".hdf":
                        data.read(filename, 0)
                    else:
                        convert_files = True
                        if os.path.exists(Path(filename).with_suffix(".hdf")):
                            message = f'\nFile\n\n {Path(filename).with_suffix(".hdf")} \n\nalready exists.'
                            message += " Previous analysis will be lost when importing."
                            message += "\n\nOverwrite existing hdf file?\n"
                            convert_files = (
                                sg.PopupYesNo(message, title="Warning !") == "Yes"
                            )
                        if convert_files:
                            if file_type == ".h5":
                                tio.hdf_data.convert(filename, None)
                                data.read(Path(filename).with_suffix(".hdf"), 0)
                            elif file_type == ".tdms":
                                data = tio.import_tdms(Path(filename), data)
                            else:
                                print(
                                    f"Importing filetype {file_type} is not implemented yet"
                                )

                    window["progbar"].UpdateBar(
                        100 * ((len(filenames) - i - 1) / len(filenames))
                    )

                    update_GUI(window, data)
                    update_plot(window, plot, data)

            window["computing"].Update("Idle")
            window["progbar"].UpdateBar(0)

        elif event in ["Save", "Save all"]:
            wait_time, data = computation_stack(
                window, plot, data, "Save graph", batch=event == "Save all"
            )

        elif event in ["Current file", "All files"]:
            data.close()
            filenames = data.open_files
            if filenames:
                window["file"].update(
                    values=filenames, set_to_index=list(filenames.index(data.filename))
                )
            else:
                window["file"].update(values=[])
            labels = data.list_labels()
            if labels:
                window["label"].update(
                    values=labels, set_to_index=list(labels.index(data.label))
                )
            else:
                window["label"].update(values=[])
            update_GUI(window, data)
            update_plot(window, plot, data)

        elif event == "Export":
            wait_time, data = computation_stack(
                window, plot, data, "Export", batch=True
            )

        elif event == "Display h5":
            os.system(
                r"C:\Users\noort\AppData\Local\HDF_Group\HDFView\3.0.0\hdfview.bat "
                + window["file"].Get()
            )

        elif event in [
            "channel_x",
            "channel_y",
            "channel_x_min",
            "channel_x_max",
            "channel_y_min",
            "channel_y_max",
            "Thread",
            "model",
            "ylog",
            "xlog",
            "color",
        ]:
            update_GUI(window, data)
            update_plot(window, plot, data)

        elif event in ["file", " selected", "label"]:
            filename = values["file"]
            label = values["label"]
            if event in ["file", " selected"]:
                data.read(filename, 0)
                labels = data.list_labels(selected_only=values[" selected"])
                index = labels.index(label) if label in labels else 0
                window["label"].update(values=labels, set_to_index=index)
            data.read(filename, values["label"])
            update_GUI(window, data)
            update_plot(window, plot, data)

        elif event in ["All", "None"]:
            data.parameters[" selected"] = 1 if event == "All" else 0
            data.save(traces=False)
            update_GUI(window, data)
            update_plot(window, plot, data)

        elif event in ["Toggle", "F1:112"]:
            data.set_parameter(
                " selected", 0 if data.get_parameter(" selected") > 0 else 1
            )
            data.save(traces=False)
            update_GUI(window, data, pars_only=True)

        elif event in ["Criteria..."]:
            criteria = GUI_select_traces(window, data)
            if criteria is not None:
                wait_time, data = computation_stack(
                    window, plot, data, "Select traces", batch=True, parameters=criteria
                )

        elif event == "Auto":
            for item in [
                "channel_x_min",
                "channel_x_max",
                "channel_y_min",
                "channel_y_max",
            ]:
                window[item].update("")
            update_plot(window, plot, data)

        elif event == "Fix":
            xrange = plot["ax"].get_xlim()
            xrange = [round(xrange[0], 2), round(xrange[1], 2)]
            window["channel_x_min"].update(xrange[0])
            window["channel_x_max"].update(xrange[1])
            yrange = plot["ax"].get_ylim()
            yrange = [round(yrange[0], 2), round(yrange[1], 2)]
            window["channel_y_min"].update(yrange[0])
            window["channel_y_max"].update(yrange[1])
            update_plot(window, plot, data)

        elif event == "Contour length":
            L = data.get_parameter("L (bp)")
            if L > 1:
                for axe in ["x", "y"]:
                    if window[f"channel_{axe}"].get() == "Z (um)":
                        window[f"channel_{axe}_min"].update(f"{-0.1*L*0.34/1000:.2f}")
                        window[f"channel_{axe}_max"].update(f"{1.1*L*0.34/1000:.2f}")
                    elif window[f"channel_{axe}"].get() == "Lfree (bp)":
                        window[f"channel_{axe}_min"].update("0")
                        window[f"channel_{axe}_max"].update(f"{1.1*L:.0f}")
            update_plot(window, plot, data)

        elif event == "Copy":
            addcopyfighandler.copyfig(plot["fig"])

        elif event == "Heatmap":
            wait_time, data = computation_stack(
                window, plot, data, batch=True, task="Heatmap", parameters=None
            )

            # GUI_update(window)
            # ax = plt.gca()
            # cs = ax.collections[0]
            # data = np.asarray(cs.get_offsets()).T
            # plot['ax'].hexbin(data[0], data[1], cmap=plt.cm.Blues, vmax=3, vmin=0, extent=ax.axis())
            # plot['canvas'].draw()

        elif event in ["Z(t)", "F(Z)"]:
            if event == "F(Z)":
                targets = ["Z (um)", "F (pN)"]
            else:
                targets = ["Time (s)", "Z (um)"]
            channels = list(data.traces.columns)
            for focus, target in zip(["channel_x", "channel_y"], targets):
                if target in channels:
                    window[focus].update(
                        values=channels, set_to_index=channels.index(target)
                    )
            for item in [
                "channel_x_min",
                "channel_x_max",
                "channel_y_min",
                "channel_y_max",
            ]:
                window[item].update("")
            update_plot(window, plot, data)

        elif event == "Y(X)":
            channels = list(data.traces.columns)
            for focus, target in zip(["channel_x", "channel_y"], ["X (um)", "Y (um)"]):
                if target in channels:
                    window[focus].update(
                        values=channels, set_to_index=channels.index(target)
                    )
            for item, value in zip(
                ["channel_x_min", "channel_x_max", "channel_y_min", "channel_y_max"],
                ["-6.7", "6.7", "-5", "5"],
            ):
                window[item].update(value)
            update_plot(window, plot, data)

        elif event == "Adjust...":
            result, batch = GUI_fit_parameters(window, data, plot)
            if result is not None:
                wait_time, data = computation_stack(
                    window,
                    plot,
                    data,
                    window["model"].get(),
                    parameters=result,
                    batch=batch,
                )

        elif event == "Delete":
            pars = list(window["params"].get()[i][0] for i in values["params"])
            wait_time, data = computation_stack(
                window, plot, data, "Delete", batch=True, parameters=pars
            )

        elif event == "Play":
            wait_time, data = computation_stack(
                window, plot, data, batch=True, task="Play"
            )

        elif event in ["Z bead", "Free DNA"]:
            params = {"L (bp)": 1.0, "A (nm)": 50.0, "S (pN)": 1000.0}
            params, batch = GUI_fix_parameters(window, data, plot, params)
            wait_time, data = computation_stack(
                window, plot, data, batch=batch, task=event, parameters=params
            )

        elif event in ["Phase jump", "Ruptures"]:
            wait_time, data = computation_stack(
                window, plot, data, batch=True, task=event
            )

        elif event == "Points...":
            criteria, batch = GUI_select_points(window, data, plot)
            if criteria is not None:
                wait_time, data = computation_stack(
                    window,
                    plot,
                    data,
                    task="Select points",
                    parameters=criteria,
                    batch=batch,
                )

        elif event in ["Break", "Escape:27"]:
            waittime, data = computation_stack(window, plot, data, stop=True)

        elif event == "Test":
            update_plot(window, plot, data)

        elif event == "selection mode":
            if values["selection mode"]:
                for axis_par in [
                    "channel_x_min",
                    "channel_x_max",
                    "channel_y_min",
                    "channel_y_max",
                ]:
                    window[axis_par].update(0)
            else:
                for axis_par in [
                    "channel_x_min",
                    "channel_x_max",
                    "channel_y_min",
                    "channel_y_max",
                ]:
                    window[axis_par].update(0)
            update_plot(window, plot, data)

        elif event == "Plot":
            df = data.read_parameter_table()
            pars = list(window["params"].get()[i][0] for i in values["params"])
            color_codes = {
                "red": "#FF0000",
                "blue": "#0000FE",
                "grey": "#bcbcbc",
                "green": "#007F00",
                "orange": "#FDA70A",
                "purple": "#800080",
            }
            colors = [color_codes[values["color"]], color_codes["grey"]]

            if len(pars) == 1:
                GUI_histogram(window, plot, df[pars[0]].dropna(), df[" selected"])
            elif len(pars) > 1:
                df = df.loc[df[" selected"] == 1]
                plt.close("all")
                sns.set_palette(sns.color_palette(colors))
                sns.pairplot(
                    df, hue="filename", vars=pars, corner=True
                )._legend.remove()
                plt.show()

        # Keyboard events
        elif event in ["F11:122", "F12:123"]:
            if data.read_next_trace(
                selected_only=values[" selected"], count_up="F12" in event
            ):
                update_GUI(window, data)
                update_plot(window, plot, data)

            # elif event in ['Left:37', 'Right:39', 'MouseWheel:Up', 'MouseWheel:Down']:
            #     if 'F11' in event:
            #         focus = 'label'
            #         event = 'Left'
            #     if 'F12' in event:
            #         focus = 'label'
            #         event = 'Right'

            #     if focus in ['file', 'label', 'channel_x', 'channel_y']:
            #         if event in ['Right', 'MouseWheel:Up']:
            #             step = 1
            #         elif event in ['Left', 'MouseWheel:Down']:
            #             step = -1
            #         else:
            #             step = 0

            #         items = list(window[focus].Widget['values'])
            #         selected = window[focus].Widget.current() + step
            #         if 0 <= selected < len(items):
            #             window[focus].update(values=items, set_to_index=selected)
            #         else:
            #             if focus in ['label']:
            #                 filenames = window['file'].Widget['values']
            #                 selected = window['file'].Widget.current() + step
            #                 if 0 <= selected < len(filenames):
            #                     window['file'].update(set_to_index=selected)
            #                     data.open(window['file'].Get(), 0)
            #                     if window[' selected'].Get():
            #                         df = data.contents.loc[(data.contents['filenames'] == filename) & (
            #                             data.contents['selected'] == 1)]
            #                     else:
            #                         df = data.contents.loc[data.contents['filenames'] == filename]

            #                     new_labels = list(df['labels'].unique())

            #                     if step == -1:
            #                         next_label = len(new_labels) - 1
            #                     else:
            #                         next_label = 0

            #                     window['label'].update(
            #                         values=new_labels, set_to_index=next_label)

            #         data.open(window['file'].Get(), window['label'].Get())
            #         update_plot(window, plot, data)

            elif focus in [
                "channel_x_min",
                "channel_x_max",
                "channel_y_min",
                "channel_y_max",
            ]:
                value = adjust_value(
                    float(window[focus].get()),
                    "Right" in event,
                    kb.is_pressed("ctrl"),
                    kb.is_pressed("shift"),
                )
                if value:
                    window[focus].update(value=f"{value:.4g}")
                update_plot(window, plot, data)

            # elif focus in ['params']:
            #     pars = window['params'].get()
            #     row = values['params'][0]
            #     value = adjust_value(pars[row][1], 'Right' in event, kb.is_pressed(
            #         'ctrl'), kb.is_pressed('shift'))
            #     if value:
            #         tio.h5_write_par(data, pars[row][0], value)
            #     update_plot(window, plot, data)
            #     window['params'].Widget.selection_set(row + 1)
            #     window['params'].Widget.focus(row + 1)  # move focus

        elif event in ["Toggle", "F3:114"]:
            print(data.parameters)

        elif event in ["Toggle", "F4:115"]:
            df = tio.h5_read_par_table(
                data.filename,
                np.asarray(window["params"].get())[:, 0],
                selected=window[" selected"].get(),
            )
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                df2 = df.describe()
            df.loc["mean"] = df2.loc["mean"]
            df.loc["std"] = df2.loc["std"]
            print("\n", df, "\n")

        elif event is not None:
            if kb.is_pressed("ctrl"):
                print("Ctrl+", {event}, "Pressed")
                pass
            elif focus in [
                "channel_x_min",
                "channel_x_max",
                "channel_y_min",
                "channel_y_max",
            ]:
                try:
                    update_plot(window, plot, data)
                except:
                    pass
            else:
                print(f"\nEvent >>>>{event}")
                print("Values >>>>", values)
                print("Focus>>>>", focus)
        # except Exception as e:
        #     exc_type, exc_obj, exc_tb = sys.exc_info()
        #     message = f'{exc_type} in {os.path.split(exc_tb.tb_frame.f_code.co_filename)[-1]}, line {exc_tb.tb_lineno}, {e}'
        #     window['computing'].Update(message)
    window.close()


if __name__ == "__main__":
    # start_files = ['data_001.h5']
    # start_files = [r'C:/Users/noort/data/171107/HMfA_data_018.h5']
    # start_files = [r'C:/Users/noort/PycharmProjects/force-spectroscopy/data_006.h5',
    #                r'C:/Users/noort/PycharmProjects/force-spectroscopy/data_007.h5',
    #                r'C:/Users/noort/PycharmProjects/force-spectroscopy/data_008.h5']
    # start_files = [r'C:\Users\noort\PycharmProjects\force-spectroscopy_v2\data_007.h5']
    # start_files.append(r'D:\tmp\data_013.h5')

    # trace_editor(start_files, '6')
    # filename = r"N:\Gert-Jan Kuijntjes N-drive\Tweezers\Waad 2nd visit - Stringent selection - finished, pull to 43pN\601_CE-chromatinNOH5_20sto43pN\data_004.h5"

    screen["plot_size"] *= 2.0

    trace_editor()
    # filename = r'\\data03\pi-vannoort\Noort\Data\Lambda DNA\data_094.h5'
    # filename = Path(rf'data\data_014.hdf')
    # trace_editor(filename, '36')
    # trace_editor(filename, '1')
