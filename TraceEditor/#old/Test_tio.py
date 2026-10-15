import tkinter as tk
from tkinter import ttk
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import pandas as pd
import numpy as np

import TraceIo as tio


class HistogramGUI:
    def __init__(self, data):
        self.data = data

        self.root = tk.Tk()
        self.root.title("Parameters histogram")

        # Create the Matplotlib figure and axis
        self.fig = Figure(figsize=(5, 4), dpi=100)
        self.ax = self.fig.add_subplot(111)

        # Create a canvas to display the Matplotlib figure
        self.canvas = FigureCanvasTkAgg(self.fig, master=self.root)
        self.canvas.draw()
        self.canvas.get_tk_widget().pack(side=tk.LEFT, fill=tk.BOTH, expand=1)

        # Create a frame for the entries and button
        self.control_frame = tk.Frame(self.root)
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

        # Create entries for #bins, min, and max
        self.bins_label = tk.Label(
            self.control_frame, text="# Bins:", anchor="w", justify="left"
        )
        self.bins_label.pack(fill=tk.X)
        self.bins_label.pack()
        self.bins_entry = tk.Entry(self.control_frame, width=18)
        self.bins_entry.pack()
        self.bins_entry.bind("<Return>", self.update_histogram)

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

        self.root.mainloop()

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

        # Clear the previous plot
        self.ax.clear()

        # Plot the histogram
        self.ax.hist(
            self.data[selected_column],
            bins=bins,
            range=(min_val, max_val),
            edgecolor="grey",
            color="Lightgrey",
        )

        # Plot the histogram for selected data
        # print(self.data)
        selected_data = self.data.loc[self.data[" selected"] == 1]
        self.ax.hist(
            selected_data[selected_column],
            bins=bins,
            range=(min_val, max_val),
            edgecolor="black",
            color="blue",
            alpha=0.5,
            label="Selected Data",
        )

        self.ax.set_xlabel(selected_column)
        self.ax.set_ylabel("Count")
        self.fig.tight_layout()

        # Redraw the canvas
        self.canvas.draw()

    def get_entry_values(self):
        bins = int(self.bins_entry.get()) if self.bins_entry.get() else 10
        min_val = (
            float(self.min_entry.get()) if self.min_entry.get() else np.min(self.data)
        )
        max_val = (
            float(self.max_entry.get()) if self.max_entry.get() else np.max(self.data)
        )
        return bins, min_val, max_val


# Example usage
if __name__ == "__main__":
    # Example DataFrame
    data = pd.DataFrame(
        {
            "A": np.random.randn(1000),
            "B": np.random.randn(1000) * 2,
            "C": np.random.randn(1000) + 1,
        }
    )

    filename = rf"data\data_010.hdf"
    label = "0"
    data = tio.hdf_data()
    data.read(filename, label)
    df = data.parameters
    df.drop("global", inplace=True)
    print(df)
    # drop colums starting with "e:"
    df = df.loc[:, ~df.columns.str.startswith("e:")]
    print(df)
    gui = HistogramGUI(df)
