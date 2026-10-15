import os
from pathlib import Path
import sys

# https://github.com/astral-sh/uv/issues/7036
if "TCL_LIBRARY" not in os.environ:
    _lib = Path(sys.base_prefix) / "lib"
    _tcl = next(_lib.glob("tcl8.*"), None)
    _tk  = next(_lib.glob("tk8.*"), None)
    if _tcl: os.environ["TCL_LIBRARY"] = str(_tcl)
    if _tk:  os.environ["TK_LIBRARY"]  = str(_tk)

import matplotlib
matplotlib.use("TkAgg")

import warnings

# Suppress the pkg_resources deprecation warning from tweezepy
warnings.filterwarnings("ignore", message=".*pkg_resources is deprecated.*")

import tkinter as tk
from tkinter import Menu
import tkinter.ttk as ttk
import tkinter.filedialog as fd
from matplotlib.figure import Figure
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import multiprocessing
from icecream import ic
from pathlib import Path
import xml.etree.ElementTree as ET
import svgwrite
import os
from icecream import ic


import numpy as np
import pandas as pd
import TraceIO as tio
import MultipProcessJob as mpj
import TraceEditorGUIs as teg
import ForceSpectroscopy as fs


class NewGUI:
    def __init__(self, root):
        self.root = root
        self.selected_only = tk.BooleanVar()
        self.x_log = tk.BooleanVar()
        self.y_log = tk.BooleanVar()
        self.x_min = tk.StringVar()
        self.x_max = tk.StringVar()
        self.y_min = tk.StringVar()
        self.y_max = tk.StringVar()
        self.model = tk.StringVar()
        self.label = tk.StringVar()
        self.color = tk.StringVar()
        self.x_axis = tk.StringVar()
        self.y_axis = tk.StringVar()
        self.parameters = None
        self.fig = Figure(figsize=(5, 4), dpi=100)
        self.ax = self.fig.add_subplot(111)
        self.previous_scale = []

        self.create_menu()
        self.create_layout()

        self.root.title("Trace editor v3")
        self.root.bind("<F12>", self.change_label)
        self.root.bind("<F11>", self.change_label)
        self.root.bind("<Control-x>", self.exit_app)
        self.root.bind("<Control-e>", self.export_file)
        self.root.bind("<F10>", self.toggle_traces)
        self.root.bind("<Escape>", lambda e: self.stop_jobs())
        self.root.bind("<F1>", self.select_traces)
        self.root.bind("<Delete>", self.delete_parameter)

        self.root.geometry("+0+0")

        # Initialize figure and axes
        self.fig = Figure(figsize=(5, 4), dpi=100)
        self.ax = self.fig.add_subplot(111)

        # Create canvas for the figure
        self.figure_canvas = FigureCanvasTkAgg(self.fig, master=self.canvas)
        self.figure_canvas.draw()
        self.figure_canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=1)

        # Call draw_figure after initialization
        self.draw_figure()

        # For testing purposes; delete this line
        # self.init_for_testing()

    ###########################################################
    # All GUI elements are created below this line
    ###########################################################

    def create_menu(self):
        self.menu_def = [
            (
                "File",
                [
                    ("Open", self.open_file),
                    ("Import", self.import_file),
                    (
                        "Close",
                        [
                            ("Current file", self.close_file),
                            ("All files", self.close_all_files),
                        ],
                    ),
                    (
                        "Save",
                        [
                            ("Export (Ctrl-e)", self.export_file),
                            # ("Save all", self.save_all_files),
                            ("Graph", self.save_graph),
                        ],
                    ),
                    ("Display h5", self.display_h5),
                    ("Exit (Ctrl-x)", self.exit_app),
                ],
            ),
            (
                "Select",
                [
                    # ("Points...", self.select_points_range),
                    (
                        "Points",
                        [
                            ("Range...", self.select_points_range),
                            ("In view", self.select_points_in_view),
                        ],
                    ),
                    (
                        "Traces",
                        [
                            ("All", self.select_all_traces),
                            ("None", self.select_none_traces),
                            ("Toggle (F10)", self.toggle_traces),
                            ("Criteria...", self.select_traces),
                        ],
                    ),
                ],
            ),
            (
                "Parameters",
                [
                    ("Fit...", self.fit_parameters),
                    ("Plot", self.plot_parameters),
                    ("Delete", self.delete_parameter),
                ],
            ),
            (
                "Process",
                [
                    (
                        "Force",
                        [
                            ("Calibrate", self.calibrate_force),
                            ("Recalculate", self.recalculate_force),
                        ],
                    ),
                    ("Free DNA", self.free_dna),
                    ("Ruptures", self.ruptures),
                ],
            ),
            (
                "Plot",
                [
                    ("Copy", self.copy_plot),
                    (
                        "Type",
                        [
                            ("Z(t)", self.z_t_plot),
                            ("F(Z)", self.f_z_plot),
                            ("Y(X)", self.y_x_plot),
                        ],
                    ),
                    (
                        "Scale",
                        [
                            ("Auto", self.auto_scale),
                            ("Fix", self.fix_scale),
                            ("Contour length", self.contour_length_scale),
                        ],
                    ),
                    # ("Play", self.play_plot),
                    ("Heatmap", self.heatmap_plot),
                ],
            ),
            (
                "Help",
                [
                    ("Break (Esc)", self.break_help),
                    ("About...", self.about_help),
                    ("Test", self.test_help),
                ],
            ),
        ]

        def add_sub_menus(parent_menu, sub_menus):
            for item in sub_menus:
                if isinstance(item[1], list):
                    sub_menu = Menu(parent_menu, tearoff=0)
                    parent_menu.add_cascade(label=item[0], menu=sub_menu)
                    add_sub_menus(sub_menu, item[1])
                else:
                    parent_menu.add_command(label=item[0], command=item[1])

        menu_bar = Menu(self.root)
        menus = {}
        for main_menu, sub_menus in self.menu_def:
            menus[main_menu] = Menu(menu_bar, tearoff=0)
            menu_bar.add_cascade(label=main_menu, menu=menus[main_menu])
            add_sub_menus(menus[main_menu], sub_menus)
        self.root.config(menu=menu_bar)

    def create_layout(self):
        # Create a frame to hold the three columns
        main_frame = tk.Frame(self.root)
        main_frame.grid(row=0, column=0, sticky="nsew")

        # Configure grid weights for resizing
        self.root.grid_rowconfigure(0, weight=1)
        self.root.grid_columnconfigure(0, weight=1)
        main_frame.grid_rowconfigure(0, weight=1)
        main_frame.grid_columnconfigure(0, weight=1)
        main_frame.grid_columnconfigure(1, weight=3)
        main_frame.grid_columnconfigure(2, weight=1)

        # Create the first column frame (model and parameters)
        parameter_column_frame = tk.Frame(main_frame)
        parameter_column_frame.grid(row=0, column=0, sticky="nsew", padx=5, pady=5)

        # Create a sub-frame for the model label and drop-down box
        model_frame = tk.Frame(parameter_column_frame)
        model_frame.grid(row=0, column=0, sticky="ew", padx=5, pady=5)

        # Create a label for the model combobox
        model_label = tk.Label(parameter_column_frame, text="Model:")
        model_label.grid(row=0, column=0, sticky="nw", padx=5, pady=(10, 0))

        # Create a combobox for the model
        fitmodels = [
            "none",
            "ExponentialDrift",
            "zMagnetForce",
            "zWLC",
            "zFiber",
            "zHMf",
            "zHMf2",
        ]
        self.model = ttk.Combobox(
            parameter_column_frame, values=fitmodels, textvariable=self.model
        )
        self.model.set("none")  # Set default value
        self.model.bind("<<ComboboxSelected>>", lambda event: self.draw_figure())
        self.model.grid(row=1, column=0, sticky="ew", padx=5, pady=(0, 5))

        # Create a label for the listbox
        parameters_label = tk.Label(parameter_column_frame, text="Parameters:")
        parameters_label.grid(row=2, column=0, sticky="nw", padx=5, pady=(5, 0))

        # Create a listbox to display data
        self.parameters = tk.Listbox(parameter_column_frame)
        self.parameters.grid(row=3, column=0, sticky="news", padx=5, pady=(0, 5))

        # Ensure the rows are configured correctly
        parameter_column_frame.grid_rowconfigure(3, weight=1)
        parameter_column_frame.grid_columnconfigure(0, weight=1)

        # Create the second column frame (filename and plot)
        plot_column_frame = tk.Frame(main_frame)
        plot_column_frame.grid(row=0, column=1, sticky="nsew", padx=5, pady=5)

        # Create a sub-frame for the filename label and drop-down box
        filename_frame = tk.Frame(plot_column_frame)
        filename_frame.grid(row=0, column=0, sticky="ew", padx=0, pady=5)

        # Create the filename drop-down box and label
        filename_label = tk.Label(filename_frame, text="Filename:")
        filename_label.grid(row=0, column=0, sticky="w", padx=0, pady=(5, 0))

        self.filename = ttk.Combobox(filename_frame, values=[], width=40)
        self.filename.grid(row=1, column=0, sticky="ew", padx=0, pady=(0, 5))
        filename_frame.grid_columnconfigure(0, weight=1)
        self.filename.bind("<<ComboboxSelected>>", lambda event: self.draw_figure())

        # Create a canvas to display the plot
        self.canvas = tk.Canvas(plot_column_frame)
        self.canvas.grid(row=1, column=0, sticky="nsew", padx=0, pady=5)
        plot_column_frame.grid_rowconfigure(1, weight=1)
        plot_column_frame.grid_columnconfigure(0, weight=1)
        self.draw_figure()

        # Create a custom progress bar using Canvas
        self.progress_canvas = tk.Canvas(
            plot_column_frame, height=10, bg="white", bd=0, highlightthickness=2
        )
        self.progress_canvas.grid(
            row=3, column=0, columnspan=5, pady=1, sticky="ew", padx=0
        )
        self.progress_bar = self.progress_canvas.create_rectangle(
            0, 0, 0, 1, fill="blue"
        )

        # Create the third column frame (label, x-axis, and y-axis)
        settings_column_frame = tk.Frame(main_frame)
        settings_column_frame.grid(row=0, column=2, sticky="nsew", padx=5, pady=5)

        # Create a sub-frame for the label drop-down box
        label_frame = tk.Frame(settings_column_frame)
        label_frame.grid(row=0, column=0, sticky="ew", padx=5, pady=5)

        # Create the label drop-down box and label
        number_label = tk.Label(
            label_frame, text="Label:            ( < : F11,  > : F12)"
        )
        number_label.grid(row=0, column=0, sticky="w", padx=5, pady=(5, 0))

        self.label = ttk.Combobox(
            label_frame, values=[1, 2, 3, 4, 5], width=10, textvariable=self.label
        )
        self.label.grid(row=1, column=0, sticky="ew", padx=5, pady=(0, 5))
        self.label.bind("<<ComboboxSelected>>", lambda event: self.draw_figure())
        label_frame.grid_columnconfigure(0, weight=1)

        # Create the 'Selected' checkbox under the 'label' combobox
        selected_checkbox = tk.Checkbutton(
            settings_column_frame,
            text="Selected only",
            variable=self.selected_only,
            command=self.draw_figure,
        )
        selected_checkbox.grid(row=1, column=0, sticky="w", padx=5, pady=(0, 5))

        # Create a sub-frame for the x-axis drop-down box
        x_axis_frame = tk.Frame(settings_column_frame)
        x_axis_frame.grid(row=2, column=0, sticky="ew", padx=5, pady=0)

        # Create the x-axis drop-down box and label
        x_axis_label = tk.Label(x_axis_frame, text="X-axis:")
        x_axis_label.grid(row=0, column=0, sticky="w", padx=5, pady=(0, 5))

        self.x_axis = ttk.Combobox(
            x_axis_frame,
            values=["Option 1", "Option 2", "Option 3"],
            width=12,
            textvariable=self.x_axis,
        )
        self.x_axis.grid(row=1, column=0, sticky="ew", padx=5, pady=(0, 5))
        x_axis_frame.grid_columnconfigure(0, weight=1)
        self.x_axis.bind("<<ComboboxSelected>>", lambda event: self.draw_figure())

        # Add minx and maxx textboxes below x-axis
        minx_label = tk.Label(x_axis_frame, text="Min:")
        minx_label.grid(row=2, column=0, sticky="w", padx=5, pady=(1, 0))

        self.x_min = tk.Entry(x_axis_frame, textvariable=self.x_min)
        self.x_min.bind("<Return>", lambda event: self.draw_figure())
        self.x_min.grid(row=3, column=0, sticky="ew", padx=5, pady=(0, 1))

        maxx_label = tk.Label(x_axis_frame, text="Max:")
        maxx_label.grid(row=4, column=0, sticky="w", padx=5, pady=(1, 0))

        self.x_max = tk.Entry(x_axis_frame, textvariable=self.x_max)
        self.x_max.bind("<Return>", lambda event: self.draw_figure())
        self.x_max.grid(row=5, column=0, sticky="ew", padx=5, pady=(0, 1))

        # Below the max_x
        log_x_checkbox = tk.Checkbutton(
            x_axis_frame,
            text="Log scale",
            variable=self.x_log,
            command=self.draw_figure,
        )
        log_x_checkbox.grid(row=6, column=0, sticky="w", padx=0, pady=(0, 5))

        # Create a sub-frame for the y-axis drop-down box
        y_axis_frame = tk.Frame(settings_column_frame)
        y_axis_frame.grid(row=3, column=0, sticky="ew", padx=5, pady=5)

        # Create the y-axis drop-down box and label
        y_axis_label = tk.Label(y_axis_frame, text="Y-axis:")
        y_axis_label.grid(row=0, column=0, sticky="w", padx=5, pady=(5, 0))

        self.y_axis = ttk.Combobox(
            y_axis_frame,
            values=["Option 1", "Option 2", "Option 3"],
            width=12,
            textvariable=self.y_axis,
        )
        self.y_axis.grid(row=1, column=0, sticky="ew", padx=5, pady=(0, 5))
        self.y_axis.bind("<<ComboboxSelected>>", lambda event: self.draw_figure())
        y_axis_frame.grid_columnconfigure(0, weight=1)

        miny_label = tk.Label(y_axis_frame, text="Min:")
        miny_label.grid(row=2, column=0, sticky="w", padx=5, pady=(0, 0))

        self.y_min = tk.Entry(y_axis_frame, textvariable=self.y_min)
        self.y_min.bind("<Return>", lambda event: self.draw_figure())
        self.y_min.grid(row=3, column=0, sticky="ew", padx=5, pady=(0, 0))

        maxy_label = tk.Label(y_axis_frame, text="Max:")
        maxy_label.grid(row=4, column=0, sticky="w", padx=5, pady=(0, 0))

        self.y_max = tk.Entry(y_axis_frame, textvariable=self.y_max)
        self.y_max.bind("<Return>", lambda event: self.draw_figure())
        self.y_max.grid(row=5, column=0, sticky="ew", padx=5, pady=(0, 0))

        # Below the max_y
        log_y_checkbox = tk.Checkbutton(
            settings_column_frame,
            text="Log scale",
            variable=self.y_log,
            command=self.draw_figure,
        )
        log_y_checkbox.grid(row=4, column=0, sticky="w", padx=5, pady=(0, 5))

        # Create a combobox below the log_scale_y checkbox
        color_label = tk.Label(settings_column_frame, text="Color:")
        color_label.grid(row=5, column=0, sticky="w", padx=10, pady=(0, 0))

        self.color = ttk.Combobox(
            settings_column_frame,
            values=["red", "blue", "green", "orange", "purple", "#1D5D86", "#33848D"],
            width=12,
            textvariable=self.color,
        )
        self.color.set("blue")  # Set default value
        self.color.grid(row=6, column=0, sticky="ew", padx=10, pady=(0, 5))
        self.color.bind("<<ComboboxSelected>>", lambda event: self.draw_figure())

    def draw_figure(self, filename=None):

        def compose_figure():
            self.ax.clear()  # Clear the existing axes
            if self.x_log.get():
                self.ax.set_xscale("log")
            if self.y_log.get():
                self.ax.set_yscale("log")

            try:
                data = tio.hdf_data(self.filename.get(), self.label.get())
            except ValueError:
                return

            old_label = self.label.get()
            self.label["values"] = data.list_labels(
                selected_only=self.selected_only.get()
            )
            if old_label not in self.label["values"]:
                self.label.set(self.label["values"][0])
            else:
                self.label.set(old_label)

            # retrieve plot data
            plot = {}
            plot["x"] = data.traces[self.x_axis.get()].values.copy()
            plot["y"] = data.traces[self.y_axis.get()].values.copy()
            try:
                plot["s"] = data.traces["Selection"]
            except KeyError:
                plot["s"] = np.ones_like(plot["x"])
            plot["xlabel"] = self.x_axis.get()
            plot["ylabel"] = self.y_axis.get()
            plot["model_name"] = self.model.get()

            # Drift corrections
            for label in ["xlabel", "ylabel"]:
                if plot[label] in ["Z (um)"] and self.model.get() != "ExponentialDrift":
                    plot[label[0]] -= data.calculate_drift()
                    plot[label[0]] -= data.get_parameter("Zbead (um)", default=0)
                if plot[label] in ["X (um)", "Y (um)"]:
                    plot[label[0]] -= data.get_parameter(
                        plot[label][0] + "0 (um)", default=0
                    )

            # plot the data
            self.ax.scatter(
                plot["x"],
                plot["y"],
                s=20,
                facecolors="none",
                edgecolors="lightgrey",
                alpha=0.5,
            )
            self.ax.scatter(
                plot["x"][plot["s"] == 1],
                plot["y"][plot["s"] == 1],
                s=20,
                facecolors="none",
                edgecolors=self.color.get(),
                alpha=0.5,
            )

            selected_parameters = None
            if plot["model_name"] not in ["none"]:
                fit_model, params, names = tio.prepare_fit(
                    plot["model_name"],
                )
                for p in params:
                    if names[p] in data.list_parameters(format="names"):
                        params[p].value = data.get_parameter(names[p])
                    if p == "z_bead":
                        params[p].value = 0

                plot["model_data"] = fit_model.eval(
                    x=data.traces[names["x"]], params=params
                )
                if names["y"] == plot["xlabel"]:
                    self.ax.plot(plot["model_data"], plot["y"], color="black")
                if names["y"] == plot["ylabel"]:
                    self.ax.plot(plot["x"], plot["model_data"], color="black")
                selected_parameters = [f" R2 {plot['model_name']}"] + list(
                    names.values()
                )
            else:
                selected_parameters = None

            # plot unwrapping for zfiber model
            if self.model.get() == "zFiber":
                l2 = data.get_parameter("L2 (bp)", default=80)
                l0 = data.get_parameter("L (bp)", default=1)
                n = data.get_parameter("n", default=0)
                if "Z (um)" in [plot["xlabel"], plot["ylabel"]]:
                    A = data.get_parameter("A (nm)", default=50)
                    S = data.get_parameter("S (pN)", default=1000)
                    f = data.traces["F (pN)"].copy().values
                    for i in range(0, int(n) + 1):
                        z = fs.zWLC(f, l0 - i * l2, A, S)
                        if plot["xlabel"] == "Z (um)":
                            self.ax.plot(z, plot["y"], color="black", alpha=0.25)
                        else:
                            self.ax.plot(plot["x"], z, color="black", alpha=0.25)
                elif "L (bp)" in [plot["xlabel"], plot["ylabel"]]:
                    smoothed_l = median_filter(data.traces["L (bp)"].copy().values, 25)
                    if plot["xlabel"] == "L (bp)":
                        self.ax.plot(smoothed_l, plot["y"], color="black")
                    else:
                        self.ax.plot(plot["x"], smoothed_l, color="black")

                    for i in range(0, int(n) + 1):
                        l = np.ones_like(plot["x"]) * (l0 - i * l2)
                        if plot["xlabel"] == "L (bp)":
                            self.ax.plot(l, plot["y"], color="black", alpha=0.25)
                        else:
                            self.ax.plot(plot["x"], l, color="black", alpha=0.25)

            self.ax.tick_params(axis="both", direction="in", top=True, right=True)
            self.ax.set_xlabel(self.x_axis.get())
            self.ax.set_ylabel(self.y_axis.get())
            self.fig.tight_layout()

            # update the listbox with the parameters
            self.parameters.delete(0, tk.END)
            for item in data.list_parameters(selection=selected_parameters):
                self.parameters.insert(tk.END, item)

            limits = [
                self.x_min.get(),
                self.x_max.get(),
                self.y_min.get(),
                self.y_max.get(),
            ]
            limits = [float(l) if l else None for l in limits]
            self.ax.set_xlim(limits[:2])
            self.ax.set_ylim(limits[2:])

            if filename:
                self.fig.savefig(filename)
            return self.fig

        # Clear the canvas
        for widget in self.canvas.winfo_children():
            widget.destroy()
        # try:
        figure = compose_figure()
        # except NameError:
        #     print("Plot not set up correctly: NameError")
        #     return

        # Create a new FigureCanvasTkAgg widget
        self.figure_canvas = FigureCanvasTkAgg(figure, master=self.canvas)
        self.figure_canvas.draw()
        self.figure_canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=1)
        self.figure_canvas.mpl_connect("button_press_event", self.on_mouse_down)
        self.figure_canvas.mpl_connect("button_release_event", self.on_mouse_up)

    def change_label(self, event):
        if event.keysym == "F11":
            # Try to decrease the label combobox value
            values = self.label["values"]
            current_index = values.index(self.label.get())
            if current_index > 0:
                self.label.set(values[current_index - 1])
                self.draw_figure()
                return
            # If label combobox cannot be increased, try to decrease the filename combobox value
            values = self.filename["values"]
            current_index = values.index(self.filename.get())
            if current_index > 0:
                self.filename.set(values[current_index - 1])
                self.label["values"] = tio.hdf_data(
                    self.filename.get(), "0"
                ).list_labels(selected_only=self.selected_only.get())
                self.label.set(self.label["values"][-1])
                self.draw_figure()
                return
        elif event.keysym == "F12":
            values = self.label["values"]
            current_index = values.index(self.label.get())
            if current_index < len(values) - 1:
                self.label.set(values[current_index + 1])
                self.draw_figure()
                return

            # If label combobox cannot be increased, try to increase the filename combobox value
            values = self.filename["values"]
            current_index = values.index(self.filename.get())
            if current_index < len(values) - 1:
                self.filename.set(values[current_index + 1])
                self.label["values"] = tio.hdf_data(
                    self.filename.get(), "0"
                ).list_labels(selected_only=self.selected_only.get())
                self.label.set(self.label["values"][0])
                self.draw_figure()
                return

    def init_for_testing(self):
        files = [
            # r"C:\tmp\vsc\TraceEditor\data\data_010.hdf",
            r"C:/tmp/vsc/TraceEditor\data\data_031.hdf",
        ]
        self.open_file(files)
        self.model.set("zMagnetForce")

    ###########################################################
    # Helper functions
    ###########################################################
    def on_mouse_down(self, event):
        # Read the coordinates of the mouse click
        if event.button == 1 and event.inaxes and event.key == "control":
            self.mouse_down = (event.xdata, event.ydata)

    def on_mouse_up(self, event):
        if event.button == 3 and event.inaxes and event.key == "control":
            # Right click with control key pressed
            if self.previous_scale:
                limits = self.previous_scale.pop()
                self.x_min.delete(0, tk.END)
                self.x_max.delete(0, tk.END)
                self.y_min.delete(0, tk.END)
                self.y_max.delete(0, tk.END)
                self.x_min.insert(0, limits[0])
                self.x_max.insert(0, limits[1])
                self.y_min.insert(0, limits[2])
                self.y_max.insert(0, limits[3])
        elif event.button == 1 and event.inaxes and event.key == "control":
            # Left click with control key pressed
            self.previous_scale.append(
                [
                    self.x_min.get(),
                    self.x_max.get(),
                    self.y_min.get(),
                    self.y_max.get(),
                ]
            )
            mouse_up = (event.xdata, event.ydata)
            x_limits = [self.mouse_down[0], mouse_up[0]]
            if x_limits[0] == x_limits[1]:
                dx = np.diff(self.ax.get_xlim()) / 5
                x_limits = [x_limits[0] - dx, x_limits[1] + dx]
            y_limits = [self.mouse_down[1], mouse_up[1]]
            if y_limits[0] == y_limits[1]:
                dy = np.diff(self.ax.get_ylim()) / 5
                y_limits = [y_limits[0] - dy, y_limits[1] + dy]

            self.x_min.delete(0, tk.END)
            self.x_max.delete(0, tk.END)
            self.y_min.delete(0, tk.END)
            self.y_max.delete(0, tk.END)
            self.x_min.insert(0, f"{np.min(x_limits):.4g}")
            self.x_max.insert(0, f"{np.max(x_limits):.4g}")
            self.y_min.insert(0, f"{np.min(y_limits):.4g}")
            self.y_max.insert(0, f"{np.max(y_limits):.4g}")
        self.draw_figure()

    def set_progress_bar(self, fraction_done):
        canvas_width = self.progress_canvas.winfo_width()
        self.progress_canvas.coords(
            self.progress_bar, 0, 0, (1 - fraction_done) * canvas_width, 10
        )
        if fraction_done == 1:
            self.draw_figure()

    def create_jobs(self, task, parameters=None):
        jobs = []
        for f in self.filename["values"]:
            data = tio.hdf_data(f, "0")
            for label in data.list_labels(selected_only=self.selected_only.get()):
                jobs.append(
                    {
                        "task": task,
                        "filename": f,
                        "label": label,
                        "parameters": parameters,
                    }
                )
        return jobs[::-1]

    def stop_jobs(self):
        self.job_process.stop()
        self.set_progress_bar(1)

    def check_queue(self):
        n_jobs = len(self.jobs)
        try:
            job_done = self.job_process.result_queue.get_nowait()
            if job_done == []:
                self.set_progress_bar(1)
            else:
                self.jobs.remove(job_done[0])
                self.set_progress_bar(1 - (len(self.jobs) / n_jobs))
        except multiprocessing.queues.Empty:
            pass
        except ValueError:
            pass
        finally:
            if not self.job_process.stop_event.is_set():
                self.progress_canvas.after(10, self.check_queue)

    def update_channels(self):
        old_x = self.x_axis.get()
        old_y = self.y_axis.get()
        channels = []
        for filename in self.filename["values"]:
            data = tio.hdf_data(filename)
            for label in data.list_labels():
                data = tio.hdf_data(filename, label)
                channels += data.list_channels()
        channels = sorted(set(channels))
        self.x_axis["values"] = channels
        self.y_axis["values"] = channels
        if old_x in channels:
            self.x_axis.set(old_x)
        elif "Time (s)" in channels:
            self.x_axis.set("Time (s)")
        else:
            self.x_axis.set(channels[0])

        if old_y in channels:
            self.y_axis.set(old_y)
        elif "Z (um)" in channels:
            self.y_axis.set("Z (um)")
        else:
            self.y_axis.set(channels[1])

    def add_text_to_svg(self, filename):
        # Check if the file exists
        if Path(filename).exists():
            # Parse the existing SVG file
            tree = ET.parse(filename)
            root = tree.getroot()

            # Define the SVG namespace
            svg_ns = "http://www.w3.org/2000/svg"
            ET.register_namespace("", svg_ns)

            # Get the width and height of the existing SVG content
            width = root.attrib.get("width", "0")
            height = root.attrib.get("height", "0")

            if width.endswith("px"):
                width = float(width[:-2])
            elif width.endswith("pt"):
                width = float(width[:-2]) * 1.333  # Convert points to pixels
            else:
                width = float(width)

            if height.endswith("px"):
                height = float(height[:-2])
            elif height.endswith("pt"):
                height = float(height[:-2]) * 1.333  # Convert points to pixels
            else:
                height = float(height)

            # Increase the width to fit the text
            new_width = width + 100  # Adjust this value as needed
            root.set("width", f"{new_width}px")

            # Adjust the viewBox to keep the plot on the left
            viewBox = root.attrib.get("viewBox", f"0 0 {width} {height}")
            viewBox_values = viewBox.split()
            viewBox_values[2] = str(new_width)
            root.set("viewBox", " ".join(viewBox_values))

            # Add text element to the right of the graph
            text = ET.Element(
                f"{{{svg_ns}}}text", x=str(width + 10), y="20", fill="black"
            )

            # Get the multi-line text
            line = os.path.basename(os.path.dirname(filename)) + "/"
            line += os.path.basename(filename).split(".")[0] + ":" + self.label.get()
            lines = [line, " "]
            lines += self.parameters.get(0, tk.END)

            # Add each line as a tspan element
            for line in lines:
                tspan = ET.Element(
                    f"{{{svg_ns}}}tspan",
                    x=str(width - 130),
                    dy="1.2em",
                    style="font-family:Arial; font-size:10px;",
                )
                tspan.text = line
                text.append(tspan)

            root.append(text)

            # Save the modified SVG file
            tree.write(filename)
        else:
            # Create a new SVG file if it doesn't exist
            dwg = svgwrite.Drawing(filename, profile="tiny")
            dwg.add(dwg.text("Hello, SVG!", insert=(10, 20), fill="black"))
            dwg.save()

    ###########################################################
    # All menu functions are implemented bewlow this line
    ###########################################################

    def open_file(self, new_files=None):
        if not new_files:
            new_files = fd.askopenfilenames(
                title="Open files",
                filetypes=[("HDF files", "*.hdf"), ("All files", "*.*")],
            )
        if new_files:
            if type(new_files) == str:
                new_files = [new_files]
            last_filename = new_files[-1]
            self.filename["values"] = sorted(
                set(list(new_files) + list(self.filename["values"]))
            )
            self.filename.set(last_filename)

            data = tio.hdf_data(new_files[-1])

            self.update_channels()
            self.label["values"] = data.list_labels()
            self.label.set(data.list_labels()[0])
        self.draw_figure()

    def import_file(self):
        new_files = fd.askopenfilenames(
            title="Open files", filetypes=[("tdms files", "*.tdms")]
        )
        for i, file in enumerate(new_files):
            self.set_progress_bar(i / len(new_files))
            data = tio.import_tdms(file)
            self.filename["values"] = sorted(
                set([str(data.filename)] + list(self.filename["values"]))
            )
            self.filename.set(data.filename)
        self.label.set(data.list_labels()[0])
        self.update_channels()
        self.set_progress_bar(1)

    def close_file(self):
        purged_files = [
            file for file in self.filename["values"] if file != self.filename.get()
        ]
        data = tio.hdf_data(self.filename.get())
        data.close()
        self.parameters.delete(0, tk.END)
        if purged_files:
            self.filename["values"] = purged_files
            self.filename.set(purged_files[0])
            data = tio.hdf_data(purged_files[0])
            self.label["values"] = data.list_labels()
            self.label.set(data.list_labels()[0])
        else:
            self.filename["values"] = []
            self.filename.set("")
            self.label["values"] = []
            self.label.set("")
            self.parameters.delete(0, tk.END)
            self.model.set("none")
        self.draw_figure()

    def close_all_files(self):
        data = tio.hdf_data(self.filename.get())
        data.close()
        self.parameters.delete(0, tk.END)
        self.filename["values"] = []
        self.filename.set("")
        self.label["values"] = []
        self.label.set("")
        self.model.set("none")
        self.draw_figure()

    def save_all_files(self):
        print("Save all files")

    def export_file(self, event=None):
        data = tio.hdf_data(self.filename.get(), label=self.label.get())
        self.set_progress_bar(0)
        filename = data.export(label=self.label.get())
        self.set_progress_bar(1)

    def save_graph(self):
        filename = fd.asksaveasfilename(
            title="Save graph",
            filetypes=[("SVG file", "*.svg"), ("JPG file", "*.jpg")],
            initialfile=Path(self.filename.get()).stem + f"_{self.label.get()}.svg",
        )
        self.draw_figure(filename)
        print(f"Data exported to {filename}")
        if filename:
            self.add_text_to_svg(filename)

    def display_h5(self):
        print("Display h5")

    def exit_app(self, event=None):
        self.root.quit()

    def select_points_range(self):
        selection_window = teg.SetSelectPointsWindow(
            self.root,
            ["Select points", self.filename.get(), self.label.get()],
            self.draw_figure,
        )

        # Bring the SelectPoints window to the front
        selection_window.lift()
        selection_window.attributes("-topmost", True)
        selection_window.after_idle(selection_window.attributes, "-topmost", False)

        # wait for the window to close
        self.root.wait_window(selection_window)
        criteria = selection_window.get_criteria()

        if len(criteria) > 1:
            self.set_progress_bar(0)
            self.jobs = self.create_jobs("select points")
            self.job_process = mpj.CalculationProcess(
                jobs=self.jobs, parameters=criteria
            )
            self.check_queue()
        elif len(criteria) == 1:
            data = tio.hdf_data(self.filename.get(), self.label.get())
            data.select_points(criteria)
            self.draw_figure()

    def select_points_in_view(self):
        criteria = [
            {
                "channel": self.x_axis.get(),
                "min": self.x_min.get(),
                "max": self.x_max.get(),
                "invert": False,
                "derivative": False,
            },
            {
                "channel": self.y_axis.get(),
                "min": self.y_min.get(),
                "max": self.y_max.get(),
                "invert": False,
                "derivative": False,
            },
        ]

        if criteria:
            data = tio.hdf_data(self.filename.get(), self.label.get())
            data.select_points(criteria)
            limits = self.previous_scale.pop()
            self.x_min.delete(0, tk.END)
            self.x_max.delete(0, tk.END)
            self.y_min.delete(0, tk.END)
            self.y_max.delete(0, tk.END)
            self.x_min.insert(0, limits[0])
            self.x_max.insert(0, limits[1])
            self.y_min.insert(0, limits[2])
            self.y_max.insert(0, limits[3])
            self.draw_figure()

    def select_all_traces(self):
        for file in self.filename["values"]:
            data = tio.hdf_data(file)
            data.parameters[" selected"] = 1
            data.save()
        self.draw_figure()

    def select_none_traces(self):
        for file in self.filename["values"]:
            data = tio.hdf_data(file)
            data.parameters[" selected"] = 0
            data.save()
        self.draw_figure()

    def toggle_traces(self, event=None):
        data = tio.hdf_data(self.filename.get(), self.label.get())
        data.set_parameter(" selected", 1 - data.get_parameter(" selected"))
        data.save(parameters=True)
        self.draw_figure()

    def select_traces(self, event=None):
        selection_window = teg.SetSelectTracesWindow(
            self.root,
            ["Select points", self.filename.get(), self.label.get()],
            self.draw_figure,
        )
        selection_window.lift()
        selection_window.attributes("-topmost", True)
        selection_window.after_idle(selection_window.attributes, "-topmost", False)

        self.root.wait_window(selection_window)
        criteria = selection_window.get_criteria()

        if criteria:
            self.set_progress_bar(0)
            self.jobs = self.create_jobs("select traces")
            self.job_process = mpj.CalculationProcess(
                jobs=self.jobs, parameters=criteria
            )
            self.check_queue()

    def fit_parameters(self):
        if self.model.get() != "none":
            fit_window = teg.SetFitParametersWindow(
                self.root,
                [self.model.get(), self.filename.get(), self.label.get()],
                self.draw_figure,
            )
            fit_window.lift()
            fit_window.attributes("-topmost", True)
            fit_window.after_idle(fit_window.attributes, "-topmost", False)
            self.root.wait_window(fit_window)

            if fit_window.button == "OK":
                self.jobs = [
                    {
                        "task": self.model.get(),
                        "filename": self.filename.get(),
                        "label": self.label.get(),
                    }
                ]
            elif fit_window.button == "Batch":
                self.jobs = self.create_jobs(self.model.get())
            else:
                self.jobs = None

            if self.jobs:
                self.set_progress_bar(0)
                self.job_process = mpj.CalculationProcess(
                    jobs=self.jobs,
                    parameters=fit_window.params,
                    fitmodels=list(self.model["values"]),
                )
                self.check_queue()
                self.draw_figure()

    def plot_parameters(self):
        data = tio.hdf_data(self.filename.get(), self.label.get())
        teg.HistogramGUI(self.root, data.parameters, color=self.color.get())

    def recalculate_force(self):

        fit_window = teg.SetFitParametersWindow(
            self.root,
            ["MagnetForce", self.filename.get(), self.label.get()],
            self.draw_figure,
        )

        fit_window.lift()
        fit_window.attributes("-topmost", True)
        fit_window.after_idle(fit_window.attributes, "-topmost", False)
        self.root.wait_window(fit_window)

        if fit_window.button == "OK":
            self.jobs = [
                {
                    "task": "zMagnetForce",
                    "filename": self.filename.get(),
                    "label": self.label.get(),
                    "parameters": fit_window.params,
                }
            ]
        elif fit_window.button == "Batch":
            self.jobs = self.create_jobs("zMagnetForce", fit_window.params)
        else:
            self.jobs = None

        if self.jobs:
            self.set_progress_bar(0)
            self.job_process = mpj.CalculationProcess(
                jobs=self.jobs,
                parameters=None,
                fitmodels=None,
            )
            self.check_queue()
            self.draw_figure()

    def calibrate_force(self, event=None):
        fit_window = teg.SetFitParametersWindow(
            self.root,
            ["power_spectrum", self.filename.get(), self.label.get()],
            self.draw_figure,
        )
        fit_window.lift()
        fit_window.attributes("-topmost", True)
        fit_window.after_idle(fit_window.attributes, "-topmost", False)
        self.root.wait_window(fit_window)

        if fit_window.button == "OK":
            self.jobs = [
                {
                    "task": "calibrate force",
                    "filename": self.filename.get(),
                    "label": self.label.get(),
                    "parameters": fit_window.params,
                }
            ]
        elif fit_window.button == "Batch":
            self.jobs = self.create_jobs("calibrate force", fit_window.params)
        else:
            self.jobs = None

        if self.jobs:
            self.set_progress_bar(0)
            self.job_process = mpj.CalculationProcess(
                jobs=self.jobs,
                parameters=None,
                fitmodels=None,
            )
            self.check_queue()
            self.draw_figure()

    def delete_parameter(self, event=None):
        parameter = self.parameters.get(tk.ACTIVE).split(" = ")[0]
        if parameter != " selected":
            data = tio.hdf_data(self.filename.get(), self.label.get())
            data.set_parameter(parameter, np.nan)
            data.save(parameters=True)
            self.draw_figure()

    def z_bead(self):
        print("Z bead")

    def free_dna(self):
        # self.jobs = self.create_jobs("Calc free DNA")
        # self.set_progress_bar(0)
        # self.job_process = mpj.CalculationProcess(jobs=self.jobs)
        # self.check_queue()
        self.set_progress_bar(0)
        data = tio.hdf_data(self.filename.get(), self.label.get())
        data.calc_free_DNA()
        data.save(parameters=False)
        self.update_channels()
        self.set_progress_bar(1)
        self.y_axis.set("L (bp)")
        self.draw_figure()

    def phase_jump(self):
        print("Phase jump")

    def ruptures(self):
        print("Ruptures")

    def copy_plot(self):
        print("Copy plot")

    def z_t_plot(self):
        self.x_axis.set("Time (s)")
        self.y_axis.set("Z (um)")
        self.contour_length_scale()
        self.draw_figure()

    def f_z_plot(self):
        self.x_axis.set("Z (um)")
        self.y_axis.set("F (pN)")
        self.contour_length_scale()
        self.draw_figure()

    def y_x_plot(self):
        self.x_axis.set("X (um)")
        self.y_axis.set("Y (um)")
        self.x_min.delete(0, tk.END)
        self.x_max.delete(0, tk.END)
        self.y_min.delete(0, tk.END)
        self.y_max.delete(0, tk.END)
        self.x_min.insert(0, "-2.5")
        self.x_max.insert(0, "2.5")
        self.y_min.insert(0, "-2.5")
        self.y_max.insert(0, "2.5")
        self.draw_figure()

    def auto_scale(self):
        self.y_min.delete(0, tk.END)
        self.y_max.delete(0, tk.END)
        self.x_min.delete(0, tk.END)
        self.x_max.delete(0, tk.END)
        self.previous_scale = []
        self.draw_figure()

    def fix_scale(self):
        if self.ax:
            x_limits = self.ax.get_xlim()
            y_limits = self.ax.get_ylim()
            self.x_min.insert(0, f"{np.min(x_limits):.3g}")
            self.x_max.insert(0, f"{np.max(x_limits):.3g}")
            self.y_min.insert(0, f"{np.min(y_limits):.3g}")
            self.y_max.insert(0, f"{np.max(y_limits):.3g}")
            self.previous_scale = []
            self.draw_figure()

    def contour_length_scale(self):
        data = tio.hdf_data(self.filename.get(), self.label.get())
        z_max = 1.2 * data.get_parameter("L (bp)") * 0.34 / 1000
        if z_max > 0:
            if self.x_axis.get() == "Z (um)":
                self.y_min.delete(0, tk.END)
                self.y_max.delete(0, tk.END)
                self.x_min.delete(0, tk.END)
                self.x_max.delete(0, tk.END)
                self.x_min.insert(0, f"{-z_max/10:.3g}")
                self.x_max.insert(0, f"{z_max:.3g}")
            if self.y_axis.get() == "Z (um)":
                self.x_min.delete(0, tk.END)
                self.x_max.delete(0, tk.END)
                self.y_min.delete(0, tk.END)
                self.y_max.delete(0, tk.END)
                self.y_min.insert(0, f"{-z_max/10:.3g}")
                self.y_max.insert(0, f"{z_max:.3g}")
            self.previous_scale = []
            self.draw_figure()

    def heatmap_plot(self):
        print("Heatmap plot")

    def break_help(self):
        print("Break help")

    def about_help(self):
        print("About help")

    def test_help(self):
        self.set_progress_bar(0.5)
        filename = fd.askopenfilename(
            title="Open files",
            filetypes=[("hdf file", "*.hdf")],
        )
        data = tio.hdf_data(filename)
        ic(data.list_channels())
        shift = data.traces["Shift (mm)"]

        for file in self.filename["values"]:
            data = tio.hdf_data(file)
            data.traces["Shift (mm)"] = shift

            data.set_parameter("L1 (mm)", 1.4, type="global")
            data.set_parameter("L2 (mm)", 0.8, type="global")
            data.set_parameter("Fmin (pN)", 0.02, type="global")
            data.set_parameter("Fmax (pN)", 85, type="global")
            data.set_parameter("alpha (-)", 0.7, type="global")
            data.set_parameter("offset (pts)", 0, type="global")

            data.traces["F (pN)"] = fs.MagnetForce(
                data.traces["Shift (mm)"],
                L1=data.parameters.loc["global", "L1 (mm)"],
                L2=data.parameters.loc["global", "L2 (mm)"],
                f0=data.parameters.loc["global", "Fmin (pN)"],
                fmax=data.parameters.loc["global", "Fmax (pN)"],
                alpha=data.parameters.loc["global", "alpha (-)"],
                dt=data.parameters.loc["global", "offset (pts)"],
            )

            # Convert shared_tracenames to list if needed and add "F (pN)"
            if isinstance(data.shared_tracenames, pd.Index):
                data.shared_tracenames = data.shared_tracenames.tolist()
            elif not isinstance(data.shared_tracenames, list):
                data.shared_tracenames = list(data.shared_tracenames)

            # Add "F (pN)" if not already present
            if "F (pN)" not in data.shared_tracenames:
                data.shared_tracenames.append("F (pN)")

            data.save()
            print(f"Data saved to {file}")


if __name__ == "__main__":
    root = tk.Tk()
    app = NewGUI(root)
    root.mainloop()
