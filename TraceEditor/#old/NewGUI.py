import tkinter as tk
from tkinter import Menu
from tkinter import ttk
import tkinter.filedialog as fd
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import numpy as np
import keyboard
import multiprocessing
from icecream import ic
from time import sleep

import TraceIo as tio
import NewGUI_Windows as window_handler
import GUI_fit as fit
import MultipProcessJob as mpj

# Function to create a dictionary of all the current values in the third column of the main_frame


def main():
    def get_plot_settings(print_values=False, prompt=None):
        plot_settings = {
            "model": model.get(),
            "filenames": filename["values"],
            "filename": filename.get(),
            "label": label.get(),
            "selected": selected.get(),
            "x_axis": x_axis.get(),
            "x_min": x_min.get(),
            "x_max": x_max.get(),
            "x_log": x_log.get(),
            "y_axis": y_axis.get(),
            "y_min": y_min.get(),
            "y_max": y_max.get(),
            "y_log": y_log.get(),
            "color": color.get(),
        }

        if prompt:
            prompt += "> "
        if print_values:
            if type(print_values) == bool:
                print(prompt)
                for name, value in plot_settings.items():
                    print(f"{name}: {value}")
            elif type(print_values) == str:
                print(f"{prompt}{print_values}: {plot_settings[print_values]}")
        return plot_settings

    # Function to handle the F12 key press

    def increase_combobox_value(event):
        # Try to increase the label combobox value
        values = label["values"]
        current_index = values.index(label.get())
        if current_index < len(values) - 1:
            label.set(values[current_index + 1])
            draw_figure()
            return

        # If label combobox cannot be increased, try to increase the filename combobox value
        values = filename["values"]
        current_index = values.index(filename.get())
        if current_index < len(values) - 1:
            filename.set(values[current_index + 1])
            label["values"] = tio.hdf_data(filename.get(), "0").list_labels()
            label.set(label["values"][0])
            draw_figure()
            return

    # Function to handle the F11 key press

    def decrease_combobox_value(event):
        # Try to decrease the label combobox value
        values = label["values"]
        current_index = values.index(label.get())
        if current_index > 0:
            label.set(values[current_index - 1])
            draw_figure()
            return

        # If label combobox cannot be increased, try to decrease the filename combobox value
        values = filename["values"]
        current_index = values.index(filename.get())
        if current_index > 0:
            filename.set(values[current_index - 1])
            label["values"] = tio.hdf_data(filename.get(), "0").list_labels()
            label.set(label["values"][-1])
            draw_figure()
            return

    # Function to draw the plot

    def compose_figure():
        global ax, fig, canvas
        if "ax" in globals():
            ax.clear()  # Clear the existing axes
        else:
            fig = Figure(figsize=(5, 4), dpi=dpi)
            ax = fig.add_subplot(111)

        if x_log.get():
            ax.set_xscale("log")
        if y_log.get():
            ax.set_yscale("log")

        try:
            if x_min.get() and x_max.get():
                ax.set_xlim(float(x_min.get()), float(x_max.get()))
            elif x_min.get():
                ax.set_xlim(float(x_min.get()), np.max(x))
            elif x_max.get():
                ax.set_xlim(np.min(x), float(x_max.get()))
            if y_min.get() and y_max.get():
                ax.set_ylim(float(y_min.get()), float(y_max.get()))
            elif y_min.get():
                ax.set_ylim(float(y_min.get()), np.max(y))
            elif y_max.get():
                ax.set_ylim(np.min(y), float(y_max.get()))
        except ValueError:
            pass

        ax.set_xlabel(x_axis.get())
        ax.set_ylabel(y_axis.get())
        fig.tight_layout()

        try:
            data = tio.hdf_data(filename.get(), label.get())
        except ValueError:
            return

        plot_names = [x_axis.get(), y_axis.get()]
        x = data.traces[plot_names[0]]
        if "Z (um)" in plot_names[0]:
            x -= data.calculate_drift() + data.get_parameter("Zbead (um)")
        if plot_names[0] in ["X (um)", "Y (um)"]:
            x -= data.get_parameter(plot_names[0][0] + "0 (um)")

        y = data.traces[plot_names[1]]
        if "Z (um)" in plot_names[1]:
            y -= data.calculate_drift() + data.get_parameter("Zbead (um)")
        if plot_names[1] in ["X (um)", "Y (um)"]:
            y -= data.get_parameter(plot_names[1][0] + "0 (um)")

        if (
            (x_axis.get() == "Time (s)")
            and (y_axis.get() == "Z (um)")
            and (model.get() == "ExponentialDrift")
        ):
            # if model.get() == "ExponentialDrift":
            print(x_axis.get() == "Time (s)", y_axis.get() == "Z (um)", model.get())
            print("ExponentialDrift")
            y += data.calculate_drift()

        try:
            s = data.traces["Selection"]
        except KeyError:
            s = np.ones_like(x)

        ax.scatter(x, y, s=20, facecolors="none", edgecolors="lightgray")

        ax.scatter(
            x[s == 1], y[s == 1], s=20, facecolors="none", edgecolors=color.get()
        )

        selected_parameters = None
        # plot model
        model_name = model.get()
        if model_name not in ["none"]:
            fit_model, params, names = tio.prepare_fit(model_name)
            if names["y"] in plot_names and names["x"] in data.traces.columns:
                for p in params:
                    if names[p] in data.list_parameters(format="names"):
                        params[p].value = data.get_parameter(names[p])
                if model_name in ["zWLC"]:
                    params["z_bead"].value = 0
                model_data = fit_model.eval(x=data.traces[names["x"]], params=params)
                if plot_names[0] == names["y"]:
                    x = model_data
                if plot_names[1] == names["y"]:
                    y = model_data
                ax.plot(x, y, color="black")

            selected_parameters = [f" R2 {model_name}"] + list(names.values())

        parameters.delete(0, tk.END)
        for item in data.list_parameters(selection=selected_parameters):
            parameters.insert(tk.END, item)

        return fig

    # Function to draw the figure on the canvas

    def draw_figure():
        # Clear the canvas
        for widget in canvas.winfo_children():
            widget.destroy()

        try:
            figure = compose_figure()
        except NameError:
            return

        # Create a new FigureCanvasTkAgg widget
        figure_canvas = FigureCanvasTkAgg(figure, master=canvas)
        figure_canvas.draw()
        figure_canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=1)

    # Function to open a pop-up window to select one or more files with the extension .hdf

    def open_file(new_files=None):
        if not new_files:
            new_files = fd.askopenfilenames(
                title="Open files",
                filetypes=[("HDF files", "*.hdf"), ("All files", "*.*")],
            )
        if new_files:
            if type(new_files) == str:
                new_files = [new_files]
            last_filename = new_files[-1]
            filename["values"] = sorted(set(list(new_files) + list(filename["values"])))
            filename.set(last_filename)

            data = tio.hdf_data(new_files[-1], "0")

            # set the contents of the combobox *-axis to the list of channels
            x_axis["values"] = data.list_channels()
            if "Time (s)" in data.list_channels():
                x_axis.set("Time (s)")
            else:
                x_axis.set(data.list_channels()[0])
            y_axis["values"] = data.list_channels()
            if "Z (um)" in data.list_channels():
                y_axis.set("Z (um)")
            else:
                y_axis.set(data.list_channels()[0])

            # set the contents of the combobox labels to the list of labels
            label["values"] = data.list_labels()
            label.set(data.list_labels()[0])
        draw_figure()

    def import_file():
        print("Import file")

    def close_file():
        print("Close file")
        filename["values"] = [f for f in filename["values"] if f != filename.get()]
        try:
            filename.set(filename["values"][0])
        except IndexError:
            filename.set("")
        draw_figure()

    def close_all_files():
        print("Close all files")
        filename["values"] = []
        filename.set("")
        draw_figure()

    def save_file():
        print("Save file")

    def save_all_files():
        print("Save all files")

    def export_file():
        print("Export file")

    def display_h5():
        print("Display h5")

    def exit_app():
        root.quit()

    def select_all_traces():
        print("Select all traces")

    def select_none_traces():
        print("Select none traces")

    def toggle_traces():
        print("Toggle traces")

    def criteria_traces():
        def print_values(values):
            print("Returned values:")
            for p, v in values.items():
                print(f"{p}: {v}")

        pars = [
            p.split(" = ")[0] for p in parameters.get(0, tk.END) if "selected" not in p
        ]
        if pars:
            criteria = window_handler.SelectTracesWindow(pars, print_values)

    def create_jobs(task):
        jobs = []
        for f in filename["values"]:
            data = tio.hdf_data(f, "0")
            for label in data.list_labels(selected_only=selected.get()):
                jobs.append((task, f, label))
        return jobs[::-1]

    def points():
        def update_selection():
            draw_figure()

        selection = fit.SetSelectPointsWindow(
            root, ["Select points", filename.get(), label.get()], update_selection
        )
        print(selection.criteria)

    def fit_parameters():
        def print_values(new_parameters):
            for par in new_parameters:
                data.set_parameter(par["parameter"], par["value"])
            data.save()
            draw_figure()

        model_name = model.get()
        if model_name == "none":
            print("No model selected")
            return
        else:
            job = create_jobs(model_name)[0]
            data = tio.hdf_data(filename.get(), label.get())
            fit_parameters = fit.SetFitParametersWindow(job, print_values)
            ic(fit_parameters.button)
            # for j in create_job(model_name):
            #     print(j)

    def plot_parameters():
        print("Plot parameters")

    def delete_parameters():
        print("Delete parameters")

    def z_bead():
        print("Z bead")

    def free_dna():
        print("Free DNA")

    def phase_jump():
        print("Phase jump")

    def ruptures():
        print("Ruptures")

    def copy_plot():
        print("Copy plot")

    def type_plot():
        print("Type plot")

    def z_t_plot():
        x_axis.set("Time (s)")
        y_axis.set("Z (um)")
        contour_length_scale()

    def f_z_plot():
        x_axis.set("Z (um)")
        y_axis.set("F (pN)")
        contour_length_scale()

    def y_x_plot():
        x_axis.set("X (um)")
        y_axis.set("Y (um)")
        auto_scale()

    def auto_scale():
        x_min.delete(0, tk.END)
        x_max.delete(0, tk.END)
        y_min.delete(0, tk.END)
        y_max.delete(0, tk.END)
        draw_figure()

    def fix_scale():
        global ax
        if ax:
            x_limits = ax.get_xlim()
            y_limits = ax.get_ylim()
            x_min.insert(0, f"{np.min(x_limits):.3g}")
            x_max.insert(0, f"{np.max(x_limits):.3g}")
            y_min.insert(0, f"{np.min(y_limits):.3g}")
            y_max.insert(0, f"{np.max(y_limits):.3g}")
            draw_figure()

    def contour_length_scale():
        data = tio.hdf_data(filename.get(), label.get())
        z_max = 1.1 * data.get_parameter("L (bp)") * 0.34 / 1000
        if z_max > 0:
            if x_axis.get() == "Z (um)":
                y_min.delete(0, tk.END)
                y_max.delete(0, tk.END)
                x_min.delete(0, tk.END)
                x_max.delete(0, tk.END)
                x_min.insert(0, "0")
                x_max.insert(0, f"{z_max:.3g}")
            if y_axis.get() == "Z (um)":
                x_min.delete(0, tk.END)
                x_max.delete(0, tk.END)
                y_min.delete(0, tk.END)
                y_max.delete(0, tk.END)
                y_min.insert(0, "0")
                y_max.insert(0, f"{z_max:.3g}")
            draw_figure()

    def play_plot():

        def stop_jobs():
            job_process.stop()
            set_progress_bar(1)

        def check_queue():
            try:
                job_done = job_process.result_queue.get_nowait()
                if job_done == []:
                    set_progress_bar(1)
                else:
                    jobs.remove(job_done[0])
                    set_progress_bar(1 - (len(jobs) / n_jobs))
                    for i, widget in enumerate([model, filename, label]):
                        widget.set(job_done[0][i])
                    draw_figure()
            except multiprocessing.queues.Empty:
                pass
            except ValueError:
                pass
            finally:
                if not job_process.stop_event.is_set():
                    progress_canvas.after(10, check_queue)

        set_progress_bar(0)
        jobs = create_jobs(model.get())
        n_jobs = len(jobs)

        job_process = mpj.CalculationProcess(jobs=jobs)
        keyboard.add_hotkey("Escape", stop_jobs)

        check_queue()

    def heatmap_plot():
        print("Heatmap plot")

    def set_progress_bar(fraction_done):
        canvas_width = progress_canvas.winfo_width()
        progress_canvas.coords(
            progress_bar, 0, 0, (1 - fraction_done) * canvas_width, 10
        )

    def break_help():
        print("Break help")
        # interrupt the calculation process
        set_progress_bar(np.random.random())

    def about_help():
        print("About help")

    def test_help():
        print("Test help")
        # use draw_plot function to draw the plot
        # get_plot_settings(print_values=True)
        # fig = compose_figure()
        # Assuming value is between 0 and 100
        # draw_figure()

    fitmodels = [
        "none",
        "ExponentialDrift",
        "MagnetForce",
        "zWLC",
        "zFiber",
        "zHMf",
        "zHMf2",
    ]

    dpi = 100

    # Create the main window
    root = tk.Tk()
    root.title("Trace editor v2")

    root.bind("<F12>", increase_combobox_value)
    root.bind("<F11>", decrease_combobox_value)

    # Create a frame to hold the three columns
    main_frame = tk.Frame(root)
    main_frame.grid(row=0, column=0, sticky="nsew")

    # Configure grid weights for resizing
    root.grid_rowconfigure(0, weight=1)
    root.grid_columnconfigure(0, weight=1)
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
    model = ttk.Combobox(parameter_column_frame, values=fitmodels)
    model.set("none")  # Set default value
    model.bind("<<ComboboxSelected>>", lambda event: draw_figure())
    model.grid(row=1, column=0, sticky="ew", padx=5, pady=(0, 5))

    # Create a label for the listbox
    parameters_label = tk.Label(parameter_column_frame, text="Parameters:")
    parameters_label.grid(row=2, column=0, sticky="nw", padx=5, pady=(5, 0))

    # Create a listbox to display data
    parameters = tk.Listbox(parameter_column_frame)
    parameters.grid(row=3, column=0, sticky="news", padx=5, pady=(0, 5))

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

    filename = ttk.Combobox(filename_frame, values=[], width=100)
    filename.grid(row=1, column=0, sticky="ew", padx=0, pady=(0, 5))
    filename_frame.grid_columnconfigure(0, weight=1)
    filename.bind("<<ComboboxSelected>>", lambda event: draw_figure())

    # Create a canvas to display the plot
    canvas = tk.Canvas(plot_column_frame)
    canvas.grid(row=1, column=0, sticky="nsew", padx=0, pady=5)
    plot_column_frame.grid_rowconfigure(1, weight=1)
    plot_column_frame.grid_columnconfigure(0, weight=1)

    # Create an empty plot initially
    empty_fig = Figure(figsize=(5, 4), dpi=100)
    empty_ax = empty_fig.add_subplot(111)
    draw_figure()

    # Create a custom progress bar using Canvas
    progress_canvas = tk.Canvas(
        plot_column_frame, height=10, bg="white", bd=0, highlightthickness=2
    )
    progress_canvas.grid(row=3, column=0, columnspan=5, pady=1, sticky="ew", padx=0)
    progress_bar = progress_canvas.create_rectangle(0, 0, 0, 1, fill="blue")

    # Create the third column frame (label, x-axis, and y-axis)
    settings_column_frame = tk.Frame(main_frame)
    settings_column_frame.grid(row=0, column=2, sticky="nsew", padx=5, pady=5)

    # Create a sub-frame for the label drop-down box
    label_frame = tk.Frame(settings_column_frame)
    label_frame.grid(row=0, column=0, sticky="ew", padx=5, pady=5)

    # Create the label drop-down box and label
    number_label = tk.Label(label_frame, text="Label:")
    number_label.grid(row=0, column=0, sticky="w", padx=5, pady=(5, 0))

    label = ttk.Combobox(label_frame, values=[1, 2, 3, 4, 5], width=10)
    label.grid(row=1, column=0, sticky="ew", padx=5, pady=(0, 5))
    label.bind("<<ComboboxSelected>>", lambda event: draw_figure())
    label_frame.grid_columnconfigure(0, weight=1)

    # Create the 'Selected' checkbox under the 'label' combobox
    selected = tk.BooleanVar()
    selected_checkbox = tk.Checkbutton(
        settings_column_frame, text="Selected only", variable=selected
    )
    selected_checkbox.grid(row=1, column=0, sticky="w", padx=5, pady=(0, 5))

    # Create a sub-frame for the x-axis drop-down box
    x_axis_frame = tk.Frame(settings_column_frame)
    x_axis_frame.grid(row=2, column=0, sticky="ew", padx=5, pady=0)

    # Create the x-axis drop-down box and label
    x_axis_label = tk.Label(x_axis_frame, text="X-axis:")
    x_axis_label.grid(row=0, column=0, sticky="w", padx=5, pady=(0, 5))

    x_axis = ttk.Combobox(
        x_axis_frame, values=["Option 1", "Option 2", "Option 3"], width=20
    )
    x_axis.grid(row=1, column=0, sticky="ew", padx=5, pady=(0, 5))
    x_axis_frame.grid_columnconfigure(0, weight=1)
    x_axis.bind("<<ComboboxSelected>>", lambda event: draw_figure())

    # Add minx and maxx textboxes below x-axis
    minx_label = tk.Label(x_axis_frame, text="Min:")
    minx_label.grid(row=2, column=0, sticky="w", padx=5, pady=(1, 0))

    x_min = tk.Entry(x_axis_frame)
    x_min.bind("<Return>", lambda event: draw_figure())
    x_min.grid(row=3, column=0, sticky="ew", padx=5, pady=(0, 1))

    maxx_label = tk.Label(x_axis_frame, text="Max:")
    maxx_label.grid(row=4, column=0, sticky="w", padx=5, pady=(1, 0))

    x_max = tk.Entry(x_axis_frame)
    x_max.bind("<Return>", lambda event: draw_figure())
    x_max.grid(row=5, column=0, sticky="ew", padx=5, pady=(0, 1))

    # Below the max_x
    x_log = tk.BooleanVar()
    log_x_checkbox = tk.Checkbutton(
        x_axis_frame, text="Log scale", variable=x_log, command=draw_figure
    )
    log_x_checkbox.grid(row=6, column=0, sticky="w", padx=0, pady=(0, 5))

    # Create a sub-frame for the y-axis drop-down box
    y_axis_frame = tk.Frame(settings_column_frame)
    y_axis_frame.grid(row=3, column=0, sticky="ew", padx=5, pady=5)

    # Create the y-axis drop-down box and label
    y_axis_label = tk.Label(y_axis_frame, text="Y-axis:")
    y_axis_label.grid(row=0, column=0, sticky="w", padx=5, pady=(5, 0))

    y_axis = ttk.Combobox(
        y_axis_frame, values=["Option 1", "Option 2", "Option 3"], width=20
    )
    y_axis.grid(row=1, column=0, sticky="ew", padx=5, pady=(0, 5))
    y_axis.bind("<<ComboboxSelected>>", lambda event: draw_figure())
    y_axis_frame.grid_columnconfigure(0, weight=1)

    miny_label = tk.Label(y_axis_frame, text="Min:")
    miny_label.grid(row=2, column=0, sticky="w", padx=5, pady=(0, 0))

    y_min = tk.Entry(y_axis_frame)
    y_min.bind("<Return>", lambda event: draw_figure())
    y_min.grid(row=3, column=0, sticky="ew", padx=5, pady=(0, 0))

    maxy_label = tk.Label(y_axis_frame, text="Max:")
    maxy_label.grid(row=4, column=0, sticky="w", padx=5, pady=(0, 0))

    y_max_var = tk.StringVar()
    y_max = tk.Entry(y_axis_frame)
    y_max.bind("<Return>", lambda event: draw_figure())
    y_max.grid(row=5, column=0, sticky="ew", padx=5, pady=(0, 0))

    # Below the max_y
    y_log = tk.BooleanVar()
    log_y_checkbox = tk.Checkbutton(
        settings_column_frame, text="Log scale", variable=y_log, command=draw_figure
    )
    log_y_checkbox.grid(row=4, column=0, sticky="w", padx=5, pady=(0, 5))

    # Create a combobox below the log_scale_y checkbox
    color_label = tk.Label(settings_column_frame, text="Color:")
    color_label.grid(row=5, column=0, sticky="w", padx=10, pady=(0, 0))

    color = ttk.Combobox(
        settings_column_frame,
        values=["red", "blue", "green", "orange", "purple"],
        width=20,
    )
    color.set("blue")  # Set default value
    color.grid(row=6, column=0, sticky="ew", padx=10, pady=(0, 5))
    color.bind("<<ComboboxSelected>>", lambda event: draw_figure())

    # Create a menubar
    menubar = Menu(root)

    # Define the menu structure
    menu_def = [
        (
            "File",
            [
                ("Open", open_file),
                ("Import", import_file),
                (
                    "Close",
                    [
                        ("Current file", close_file),
                        ("All files", close_all_files),
                    ],
                ),
                (
                    "Save",
                    [
                        ("Save all", save_all_files),
                        ("Export", export_file),
                    ],
                ),
                ("Display h5", display_h5),
                ("Exit", exit_app),
            ],
        ),
        (
            "Select",
            [
                ("Points...", points),
                (
                    "Traces",
                    [
                        ("All", select_all_traces),
                        ("None", select_none_traces),
                        ("Toggle", toggle_traces),
                        ("Criteria...", criteria_traces),
                    ],
                ),
            ],
        ),
        (
            "Parameters",
            [
                ("Fit...", fit_parameters),
                ("Plot", plot_parameters),
                ("Delete", delete_parameters),
            ],
        ),
        (
            "Process",
            [
                ("Z bead", z_bead),
                ("Free DNA", free_dna),
                ("Phase jump", phase_jump),
                ("Ruptures", ruptures),
            ],
        ),
        (
            "Plot",
            [
                ("Copy", copy_plot),
                (
                    "Type",
                    [
                        ("Z(t)", z_t_plot),
                        ("F(Z)", f_z_plot),
                        ("Y(X)", y_x_plot),
                    ],
                ),
                (
                    "Scale",
                    [
                        ("Auto", auto_scale),
                        ("Fix", fix_scale),
                        ("Contour length", contour_length_scale),
                    ],
                ),
                ("Play", play_plot),
                ("Heatmap", heatmap_plot),
            ],
        ),
        (
            "Help",
            [
                ("Break", break_help),
                ("About...", about_help),
                ("Test", test_help),
            ],
        ),
    ]

    # Function to create the menu
    def create_menu(menu_def, menubar):
        shortcuts = {"Open": "Ctrl+O", "Exit": "Ctrl+X"}
        for menu_label, menu_items in menu_def:
            submenu = Menu(menubar, tearoff=0)
            for item in menu_items:
                if isinstance(item, tuple) and isinstance(item[1], list):
                    subsubmenu = Menu(submenu, tearoff=0)
                    for subitem_label, subitem_command in item[1]:
                        subsubmenu.add_command(
                            label=subitem_label, command=subitem_command
                        )
                    submenu.add_cascade(label=item[0], menu=subsubmenu)
                else:
                    submenu.add_command(
                        label=item[0],
                        command=item[1],
                        accelerator=shortcuts.get(item[0]),
                    )
            menubar.add_cascade(label=menu_label, menu=submenu)

        root.bind_all("<Control-o>", lambda event: open_file())
        root.bind_all("<Control-x>", lambda event: exit_app())
        # root.bind_all("<Escape>", lambda event: break_help())

    # Create the menu

    create_menu(menu_def, menubar)

    root.config(menu=menubar)

    # set defaults
    open_file([r"C:\tmp\vsc\TraceEditor\data\data_014.hdf"])
    # get_plot_settings(print_values="filenames", prompt='1')

    filename["values"] = [
        r"C:\tmp\vsc\TraceEditor\data\data_014.hdf",
        r"C:/tmp/vsc/TraceEditor/data/data_031.hdf",
    ]
    filename["values"] = [r"C:/tmp/vsc/TraceEditor/data/data_014.hdf"]
    # get_plot_settings(print_values="filenames", prompt='2')

    # Start the Tkinter event loop
    root.mainloop()


if __name__ == "__main__":
    main()
