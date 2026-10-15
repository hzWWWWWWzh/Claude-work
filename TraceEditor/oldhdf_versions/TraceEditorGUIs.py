from tkinter import ttk
import tkinter as tk
from icecream import ic
import numpy as np
import tkinter as tk
from tkinter import ttk
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import pandas as pd

import TraceIO as tio


class SetFitParametersWindow(tk.Toplevel):
    def __init__(self, parent, job, update_func):
        super().__init__(parent)
        self.parent = parent
        self.job = job
        self.update_func = update_func
        self.title("Set Fit Parameters ...")
        self.button = None
        self.data = tio.hdf_data(job[1], job[2])

        self.model, self.params, self.names = tio.prepare_fit(job[0])
        label = job[2]
        update_func()

        for p in self.params:
            if self.names[p] not in self.data.parameters.columns:
                self.params[p].type = "Global"
                self.params[p].vary = False
            elif np.isfinite(self.data.parameters.loc[label, self.names[p]]):
                # parameter is in data.parameters
                self.params[p].value = self.data.parameters.loc[label, self.names[p]]
                self.params[p].type = "Local"
                self.params[p].vary = False
                try:
                    self.params[p].error = self.data.parameters.loc[
                        label, "e:" + self.names[p]
                    ]
                    if self.params[p].error > 0:
                        # parameter has an error => was fitted
                        self.params[p].type = "Fit"
                        self.params[p].vary = True
                        self.params[p].value = float(
                            tio.round_significance(
                                self.params[p].value, self.params[p].error
                            ).split(" ± ")[0]
                        )
                except KeyError:
                    pass
            elif np.isfinite(self.data.parameters.loc["global", self.names[p]]):
                # data is in global parameters
                self.params[p].value = self.data.parameters.loc["global", self.names[p]]
                self.params[p].type = "Global"
                self.params[p].vary = False

        self.old_data_params = self.data.parameters.copy()

        # Create header labels
        width = 15
        tk.Label(self, text="Channel", width=20, anchor="w").grid(
            row=0, column=0, padx=1, pady=1, sticky="w"
        )
        tk.Label(self, text="Value", width=width, anchor="w").grid(
            row=0, column=1, padx=1, pady=1, sticky="w"
        )
        tk.Label(self, text="Type", width=width, anchor="w").grid(
            row=0, column=2, padx=1, pady=1, sticky="w"
        )
        tk.Label(self, text="Min", width=width, anchor="w").grid(
            row=0, column=3, padx=1, pady=1, sticky="w"
        )
        tk.Label(self, text="Max", anchor="w", width=width).grid(
            row=0, column=4, padx=1, pady=1, sticky="w"
        )

        self.entries = []
        for i, p in enumerate(self.params):
            tk.Label(self, text=self.names[p], width=20).grid(
                row=i + 1, column=0, padx=1, pady=1, sticky="we"
            )

            value_entry = tk.Entry(self, width=width)
            value_entry.insert(0, self.params[p].value)
            value_entry.grid(row=i + 1, column=1, padx=1, pady=1, sticky="w")
            value_entry.bind(
                "<Up>",
                lambda event, entry=value_entry: tio.adjust_value(
                    entry, 1, event, self.apply_model
                ),
            )

            value_entry.bind(
                "<Down>",
                lambda event, entry=value_entry: tio.adjust_value(
                    entry, -1, event, self.apply_model
                ),
            )
            value_entry.bind("<Return>", lambda event: self.apply_model())

            type_combo = ttk.Combobox(
                self, values=["Global", "Local", "Fit"], width=width - 4
            )
            type_combo.set(self.params[p].type)
            type_combo.grid(row=i + 1, column=2, padx=8, pady=2, sticky="w")
            type_combo.bind("<<ComboboxSelected>>", lambda event: self.apply_model())

            min_entry = tk.Entry(self, width=width)
            min_entry.insert(
                0, self.params[p].min if np.isfinite(self.params[p].min) else ""
            )
            min_entry.grid(row=i + 1, column=3, padx=1, pady=1, sticky="w")
            min_entry.bind(
                "<Up>",
                lambda event, entry=min_entry: tio.adjust_value(
                    entry, 1, event, self.apply_model
                ),
            )

            min_entry.bind(
                "<Down>",
                lambda event, entry=min_entry: tio.adjust_value(
                    entry, -1, event, self.apply_model
                ),
            )
            min_entry.bind("<Return>", lambda event: self.apply_model())

            max_entry = tk.Entry(self, width=width)
            max_entry.insert(
                0, self.params[p].max if np.isfinite(self.params[p].max) else ""
            )
            max_entry.grid(row=i + 1, column=4, padx=1, pady=1, sticky="w")
            max_entry.bind(
                "<Up>",
                lambda event, entry=max_entry: tio.adjust_value(
                    entry, 1, self.apply_model
                ),
            )
            max_entry.bind(
                "<Down>",
                lambda event, entry=max_entry: tio.adjust_value(
                    entry, -1, event, self.apply_model
                ),
            )
            max_entry.bind("<Return>", lambda event: self.apply_model())
            self.entries.append((value_entry, type_combo, min_entry, max_entry))

        # Create a frame for the buttons
        button_frame = tk.Frame(self)
        button_frame.grid(row=i + 2, column=0, columnspan=5, pady=5)

        self.cancel_button = tk.Button(
            button_frame, text="Cancel", command=self.cancel, width=20
        )
        self.cancel_button.grid(row=0, column=0, padx=5, pady=1, sticky="ew")

        self.batch_button = tk.Button(
            button_frame, text="Batch", command=self.batch, width=20
        )
        self.batch_button.grid(row=0, column=1, padx=5, pady=1, sticky="ew")

        self.ok_button = tk.Button(button_frame, text="OK", command=self.ok, width=20)
        self.ok_button.grid(row=0, column=2, padx=5, pady=1, sticky="ew")

        # Configure the columns to expand equally
        button_frame.grid_columnconfigure(0, weight=1)
        button_frame.grid_columnconfigure(1, weight=1)
        button_frame.grid_columnconfigure(2, weight=1)
        button_frame.grid_columnconfigure(3, weight=1)

        # Update idle tasks to ensure geometry is calculated
        self.update_idletasks()

        # Schedule the geometry setting to run after a short delay
        self.after(100, self.set_window_position)

    def set_window_position(self):
        # Get window dimensions
        window_width = self.winfo_width()
        window_height = self.winfo_height()

        # Calculate screen dimensions
        screen_width = self.winfo_screenwidth()
        screen_height = self.winfo_screenheight()

        # Set window position to bottom-right corner
        x_offset = screen_width - window_width - 100
        y_offset = screen_height - window_height - 100
        self.geometry(f"{window_width}x{window_height}+{x_offset}+{y_offset}")

    def update_params(self):
        for i, p in enumerate(self.params):
            self.params[p].value = float(self.entries[i][0].get())
            type = self.entries[i][1].get()
            try:
                self.params[p].min = float(self.entries[i][2].get())
            except ValueError:
                self.params[p].min = -np.inf
            try:
                self.params[p].max = float(self.entries[i][3].get())
            except ValueError:
                self.params[p].max = np.inf

            self.params[p].error = None
            if type == "Fit":
                self.params[p].vary = True
            elif type == "Local":
                self.params[p].vary = False
                self.params[p].set(brute_step=0.1)
            else:
                self.params[p].vary = False
                self.params[p].set(brute_step=1.0)

    def apply_model(self):
        self.update_params()
        for p in self.params:
            # print(
            #     p,
            #     self.params[p].value,
            #     self.params[p].type,
            #     self.params[p].min,
            #     self.params[p].max,
            # )
            self.data.set_parameter(self.names[p], self.params[p].value, error=None)

        if np.any([p.vary for p in self.params.values()]):
            result = self.data.fit(self.job[0], self.params)
            if result.success:
                for i, p in enumerate(self.params):
                    if self.params[p].vary:
                        self.data.set_parameter(
                            self.names[p],
                            result.params[p].value,
                            result.params[p].stderr,
                        )
                    self.entries[i][0].delete(0, tk.END)
                    self.entries[i][0].insert(
                        0,
                        tio.round_significance(
                            result.params[p].value, result.params[p].stderr
                        ).split(" ± ")[0],
                    )

        self.data.save(parameters=True)
        self.update_func()

    def cancel(self):
        self.data.parameters = self.old_data_params
        self.data.save()
        self.button = "Cancel"
        self.update_func()
        self.destroy()  # Close the window and terminate the program

    def batch(self):
        self.update_params()
        self.button = "Batch"
        self.destroy()  # Close the window and terminate the program
        return

    def ok(self):
        self.update_params()
        self.button = "OK"
        self.destroy()  # Close the window and terminate the program
        return

    def run(self):
        self.mainloop()


class SetSelectPointsWindow(tk.Toplevel):
    def __init__(self, parent, job, update_func):
        super().__init__(parent)
        self.parent = parent
        self.job = job
        self.update_func = update_func
        self.title("Select points ...")
        self.button = None
        self.names = []
        self.criteria = []
        self.data = tio.hdf_data(job[1], job[2])

        try:
            self.old_selection = self.data.traces["Selection"]
        except KeyError:
            self.old_selection = np.ones_like(
                self.data.traces[self.data.list_channels()[0]]
            )
            self.data.traces["Selection"] = self.old_selection
            self.data.save()

        # Create header labels
        width = 15
        tk.Label(self, text="Channel", width=20, anchor="w").grid(
            row=0, column=0, padx=1, pady=1, sticky="w"
        )
        tk.Label(self, text="Min", width=width, anchor="w").grid(
            row=0, column=1, padx=1, pady=1, sticky="w"
        )
        tk.Label(self, text="Max", anchor="w", width=width).grid(
            row=0, column=2, padx=1, pady=1, sticky="w"
        )
        tk.Label(self, text="Derivative", width=width, anchor="w").grid(
            row=0, column=3, padx=1, pady=1, sticky="w"
        )
        tk.Label(self, text="Invert", width=width, anchor="w").grid(
            row=0, column=4, padx=1, pady=1, sticky="w"
        )

        self.entries = []
        for i, channel in enumerate(self.data.list_channels()):
            self.names.append(channel)
            tk.Label(self, text=channel, width=20).grid(
                row=i + 1, column=0, padx=1, pady=1, sticky="we"
            )

            min_entry = tk.Entry(self, width=width)
            min_entry.grid(row=i + 1, column=1, padx=1, pady=1, sticky="w")
            min_entry.bind(
                "<Up>",
                lambda event, entry=min_entry: tio.adjust_value(
                    entry, 1, event, self.apply_criteria
                ),
            )
            min_entry.bind(
                "<Down>",
                lambda event, entry=min_entry: tio.adjust_value(
                    entry, -1, event, self.apply_criteria
                ),
            )
            min_entry.bind("<Return>", lambda event: self.apply_criteria())

            max_entry = tk.Entry(self, width=width)
            max_entry.grid(row=i + 1, column=2, padx=1, pady=1, sticky="w")
            max_entry.bind(
                "<Up>",
                lambda event, entry=max_entry: tio.adjust_value(
                    entry, 1, event, self.apply_criteria
                ),
            )
            max_entry.bind(
                "<Down>",
                lambda event, entry=max_entry: tio.adjust_value(
                    entry, -1, event, self.apply_criteria
                ),
            )
            max_entry.bind("<Return>", lambda event: self.apply_criteria())

            derivative_combo = ttk.Combobox(
                self, values=["False", "True"], width=width // 2
            )
            derivative_combo.set("True" if channel == "Shift (mm)" else "False")
            derivative_combo.grid(row=i + 1, column=3, padx=8, pady=2, sticky="w")
            derivative_combo.bind(
                "<<ComboboxSelected>>", lambda event: self.apply_criteria()
            )

            invert_combo = ttk.Combobox(
                self, values=["False", "True"], width=width // 2
            )
            invert_combo.set("False")
            invert_combo.grid(row=i + 1, column=4, padx=8, pady=2, sticky="w")
            invert_combo.bind(
                "<<ComboboxSelected>>", lambda event: self.apply_criteria()
            )

            self.entries.append((min_entry, max_entry, derivative_combo, invert_combo))

        # Create a frame for the buttons
        button_frame = tk.Frame(self)
        button_frame.grid(row=i + 2, column=0, columnspan=5, pady=5)

        self.cancel_button = tk.Button(
            button_frame, text="Cancel", command=self.cancel, width=20
        )
        self.cancel_button.grid(row=0, column=0, padx=5, pady=1, sticky="ew")

        self.batch_button = tk.Button(
            button_frame, text="Batch", command=self.batch, width=20
        )
        self.batch_button.grid(row=0, column=1, padx=5, pady=1, sticky="ew")

        self.ok_button = tk.Button(button_frame, text="OK", command=self.ok, width=20)
        self.ok_button.grid(row=0, column=2, padx=5, pady=1, sticky="ew")

        # Configure the columns to expand equally
        button_frame.grid_columnconfigure(0, weight=1)
        button_frame.grid_columnconfigure(1, weight=1)
        button_frame.grid_columnconfigure(2, weight=1)
        button_frame.grid_columnconfigure(3, weight=1)

        # Update idle tasks to ensure geometry is calculated
        self.update_idletasks()

        # Schedule the geometry setting to run after a short delay
        self.after(100, self.set_window_position)

    def set_window_position(self):
        # Get window dimensions
        window_width = self.winfo_width()
        window_height = self.winfo_height()

        # Calculate screen dimensions
        screen_width = self.winfo_screenwidth()
        screen_height = self.winfo_screenheight()

        # Set window position to bottom-right corner
        x_offset = screen_width - window_width - 100
        y_offset = screen_height - window_height - 100
        self.geometry(f"{window_width}x{window_height}+{x_offset}+{y_offset}")

    def apply_criteria(self):
        self.criteria = []
        for i, name in enumerate(self.names):
            self.criteria.append(
                {
                    "channel": name,
                    "min": self.entries[i][0].get(),
                    "max": self.entries[i][1].get(),
                    "derivative": self.entries[i][2].get() == "True",
                    "invert": self.entries[i][3].get() == "True",
                }
            )

        self.data.select_points(self.criteria)
        self.update_func()

    def cancel(self):
        self.data.traces["Selection"] = self.old_selection
        self.data.save()
        self.button = "Cancel"
        self.update_func()
        self.criteria = []
        self.destroy()  # Close the window and terminate the program

    def batch(self):
        self.criteria = []
        for i, name in enumerate(self.names):
            self.criteria.append(
                {
                    "channel": name,
                    "min": self.entries[i][0].get(),
                    "max": self.entries[i][1].get(),
                    "derivative": self.entries[i][2].get() == "True",
                    "invert": self.entries[i][3].get() == "True",
                }
            )
        self.criteria = [c for c in self.criteria if c["min"] != "" or c["max"] != ""]
        self.button = "Batch"
        self.destroy()  # Close the window and terminate the program

    def ok(self):
        self.button = "OK"
        self.criteria = []
        self.destroy()  # Close the window and terminate the program

    def get_criteria(self):
        return self.criteria

    def run(self):
        self.mainloop()


class SetSelectTracesWindow(tk.Toplevel):
    def __init__(self, parent, job, update_func):
        super().__init__(parent)
        self.parent = parent
        self.job = job
        self.update_func = update_func
        self.title("Select points ...")
        self.button = None
        self.names = []
        self.criteria = []
        self.data = tio.hdf_data(job[1], job[2])

        # Create header labels
        width = 15
        tk.Label(self, text="Parameter", width=20, anchor="w").grid(
            row=0, column=0, padx=1, pady=1, sticky="w"
        )
        tk.Label(self, text="Min", width=width, anchor="w").grid(
            row=0, column=1, padx=1, pady=1, sticky="w"
        )
        tk.Label(self, text="Max", anchor="w", width=width).grid(
            row=0, column=2, padx=1, pady=1, sticky="w"
        )
        tk.Label(self, text="Invert", width=width, anchor="w").grid(
            row=0, column=3, padx=1, pady=1, sticky="w"
        )

        self.entries = []
        for i, channel in enumerate(self.data.list_parameters()):
            self.names.append(channel.split(" = ")[0])
            tk.Label(self, text=self.names[-1], width=20).grid(
                row=i + 1, column=0, padx=1, pady=1, sticky="we"
            )

            min_entry = tk.Entry(self, width=width)
            min_entry.grid(row=i + 1, column=1, padx=1, pady=1, sticky="w")

            max_entry = tk.Entry(self, width=width)
            max_entry.grid(row=i + 1, column=2, padx=1, pady=1, sticky="w")

            invert_combo = ttk.Combobox(
                self, values=["False", "True"], width=width // 2
            )
            invert_combo.set("False")
            invert_combo.grid(row=i + 1, column=3, padx=8, pady=2, sticky="w")

            self.entries.append((min_entry, max_entry, invert_combo))

        # Create a frame for the buttons
        button_frame = tk.Frame(self)
        button_frame.grid(row=i + 2, column=0, columnspan=5, pady=5)

        self.cancel_button = tk.Button(
            button_frame, text="Cancel", command=self.cancel, width=20
        )
        self.cancel_button.grid(row=0, column=0, padx=5, pady=1, sticky="ew")

        self.ok_button = tk.Button(button_frame, text="OK", command=self.ok, width=20)
        self.ok_button.grid(row=0, column=1, padx=5, pady=1, sticky="ew")

        # Configure the columns to expand equally
        button_frame.grid_columnconfigure(0, weight=1)
        button_frame.grid_columnconfigure(1, weight=1)

        # Update idle tasks to ensure geometry is calculated
        self.update_idletasks()

        # Schedule the geometry setting to run after a short delay
        self.after(100, self.set_window_position)

    def set_window_position(self):
        # Get window dimensions
        window_width = self.winfo_width()
        window_height = self.winfo_height()

        # Calculate screen dimensions
        screen_width = self.winfo_screenwidth()
        screen_height = self.winfo_screenheight()

        # Set window position to bottom-right corner
        x_offset = screen_width - window_width - 100
        y_offset = screen_height - window_height - 100
        self.geometry(f"{window_width}x{window_height}+{x_offset}+{y_offset}")

    def cancel(self):
        self.button = "Cancel"
        self.criteria = []
        self.destroy()  # Close the window and terminate the program

    def ok(self):
        self.button = "OK"
        self.criteria = []
        for i, name in enumerate(self.names):
            self.criteria.append(
                {
                    "Parameter": name,
                    "min": self.entries[i][0].get(),
                    "max": self.entries[i][1].get(),
                    "invert": self.entries[i][2].get() == "True",
                }
            )
        self.criteria = [c for c in self.criteria if c["min"] != "" or c["max"] != ""]
        self.destroy()  # Close the window and terminate the program

    def get_criteria(self):
        return self.criteria

    def run(self):
        self.mainloop()


class HistogramGUI(tk.Toplevel):
    def __init__(self, parent, parameters_df, color="blue"):
        super().__init__(parent)
        self.data = parameters_df
        self.data.drop("global", inplace=True)
        self.data = self.data.loc[:, ~self.data.columns.str.startswith("e:")]
        self.title("Parameters histogram")
        self.color = color

        # Create the Matplotlib figure and axis
        self.fig = Figure(figsize=(5, 4), dpi=100)
        self.ax = self.fig.add_subplot(111)

        # Create a canvas to display the Matplotlib figure
        self.canvas = FigureCanvasTkAgg(self.fig, master=self)
        self.canvas.draw()
        self.canvas.get_tk_widget().pack(side=tk.LEFT, fill=tk.BOTH, expand=1)

        # Create a frame for the entries and button
        self.control_frame = tk.Frame(self)
        self.control_frame.pack(side=tk.RIGHT, fill=tk.Y, padx=10, pady=10)

        # Create a drop-down box for selecting the column
        # self.column_label = tk.Label(self.control_frame, text="Parameter")
        self.column_label = tk.Label(
            self.control_frame, text="Parameter", anchor="w", justify="left"
        )
        self.column_label.pack(fill=tk.X)
        self.column_label.pack()
        self.column_var = tk.StringVar()
        self.column_dropdown = ttk.Combobox(
            self.control_frame, textvariable=self.column_var, width=15
        )
        # only use columns with more thand 2 unique values
        used_parameters = list(self.data.loc[:, self.data.nunique() > 2].columns)
        self.column_dropdown["values"] = used_parameters
        self.column_dropdown.pack()
        self.column_dropdown.bind("<<ComboboxSelected>>", self.update_histogram)
        self.column_dropdown.current(0)

        # # Create entries for #bins, min, and max
        self.bins_label = tk.Label(
            self.control_frame, text="# Bins:", anchor="w", justify="left"
        )
        self.bins_label.pack(fill=tk.X)
        self.bins_label.pack()
        self.bins_entry = tk.Entry(self.control_frame, width=18)
        self.bins_entry.pack()
        self.bins_entry.bind("<Return>", self.update_histogram)

        # Create entries for #bins, min, and max
        self.range = tk.Label(
            self.control_frame, text="Range:", anchor="w", justify="left"
        )
        self.range.pack(fill=tk.X)
        self.range.pack()
        self.range_entry = tk.Entry(self.control_frame, width=18)
        self.range_entry.pack()
        self.range_entry.bind("<Return>", self.update_histogram)

        self.min_label = tk.Label(
            self.control_frame, text="Min:", anchor="w", justify="left"
        )
        self.min_label.pack(fill=tk.X)
        self.min_label.pack()
        self.min_entry = tk.Entry(self.control_frame, width=18)
        self.min_entry.pack()
        self.min_entry.bind("<Return>", self.update_histogram)

        self.max_label = tk.Label(
            self.control_frame, text="Max:", anchor="w", justify="left"
        )
        self.max_label.pack(fill=tk.X)
        self.max_label.pack()
        self.max_entry = tk.Entry(self.control_frame, width=18)
        self.max_entry.pack()
        self.max_entry.bind("<Return>", self.update_histogram)

        # Create a button to update the histogram
        # self.update_button = tk.Button(
        #     self.control_frame, text="Update Histogram", command=self.update_histogram
        # )
        # self.update_button.pack(pady=10)

        # Initial plot
        self.update_histogram()

        self.mainloop()

    def update_histogram(self, event=None):
        # Get the selected column
        selected_column = self.column_var.get()
        if not selected_column:
            return

        # Get the values from the entries
        bins = int(self.bins_entry.get()) if self.bins_entry.get() else 10
        min_val = (
            float(self.min_entry.get())
            if self.min_entry.get()
            else self.data[selected_column].min()
        )
        max_val = (
            float(self.max_entry.get())
            if self.max_entry.get()
            else self.data[selected_column].max()
        )

        # Clear the entire figure
        self.fig.clf()
        self.ax = self.fig.add_subplot(111)

        # Plot the histogram
        histogram, edges = np.histogram(
            self.data[selected_column], bins=bins, range=(min_val, max_val)
        )
        bins = 0.5 * (edges[:-1] + edges[1:])
        binwidth = bins[1] - bins[0]
        divider = np.sum(histogram) * binwidth

        range_max = self.range_entry.get()
        if range_max == "":
            range_max = 1.3 * np.max(histogram) / divider

        else:
            range_max = float(range_max)

        self.ax.bar(
            bins,
            histogram / divider,
            color="lightgrey",
            edgecolor="grey",
            label="All Data",
            width=binwidth,
        )
        self.ax.set_ylabel("pdf")
        self.ax.set_ylim(0, range_max)

        # Create a second y-axis for the count if it doesn't exist
        self.ax2 = self.ax.twinx()

        # Set the limits for the secondary y-axis based on the primary y-axis limits
        self.ax2.set_ylim(0, range_max * divider)
        self.ax2.set_ylabel("Count")

        # Plot the histogram for selected data
        selected_data = self.data.loc[self.data[" selected"] == 1]
        histogram, edges = np.histogram(
            selected_data[selected_column], bins=edges, range=(min_val, max_val)
        )

        self.ax.bar(
            bins,
            histogram / divider,
            color=self.color,
            edgecolor="black",
            alpha=0.5,
            width=binwidth,
        )

        self.ax.set_xlabel(selected_column)
        self.fig.tight_layout()

        # Redraw the canvas
        self.canvas.draw()

    def get_entry_values(self):
        bins = int(self.range_entry.get()) if self.range_entry.get() else 10
        min_val = (
            float(self.min_entry.get()) if self.min_entry.get() else np.min(self.data)
        )
        max_val = (
            float(self.max_entry.get()) if self.max_entry.get() else np.max(self.data)
        )
        return bins, min_val, max_val


def update_func(values=None):
    if values is not None:
        print("update_func called with values:")
        for v in values:
            print(v)
    else:
        print("update_func called without values")

