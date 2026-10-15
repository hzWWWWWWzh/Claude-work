import tkinter as tk
from tkinter import ttk
from tkinter.simpledialog import askinteger
import multiprocessing

import TraceIo as tio


class SelectTracesWindow:
    def __init__(self, parameters, update_func):
        self.win = tk.Toplevel()
        self.win.title("Selection criteria")
        self.update_func = update_func
        self.parameters = parameters

        # Create header labels
        tk.Label(self.win, text="Parameter", width=20, anchor="w").grid(
            row=0, column=0, padx=5, sticky="we"
        )
        tk.Label(self.win, text="Min", width=13, anchor="w").grid(
            row=0, column=1, padx=5, sticky="w"
        )
        tk.Label(self.win, text="Max", width=14, anchor="w").grid(
            row=0, column=2, padx=5, sticky="w"
        )
        tk.Label(self.win, text="Invert", width=5, anchor="w").grid(
            row=0, column=3, padx=5, sticky="w"
        )

        self.entries = []
        for i, param in enumerate(parameters):
            tk.Label(self.win, text=param, width=20).grid(
                row=i + 1, column=0, padx=5, pady=4, sticky="we"
            )

            min_entry = tk.Entry(self.win, width=13)
            min_entry.grid(row=i + 1, column=1, padx=5, pady=4, sticky="w")

            max_entry = tk.Entry(self.win, width=14)
            max_entry.grid(row=i + 1, column=2, padx=5, pady=4, sticky="w")

            invert_var = tk.BooleanVar()
            invert_check = tk.Checkbutton(self.win, variable=invert_var)
            invert_check.grid(row=i + 1, column=3, padx=5, pady=2, sticky="w")

            self.entries.append((min_entry, max_entry, invert_var))

        # Create a frame for the buttons
        button_frame = tk.Frame(self.win)
        button_frame.grid(row=len(parameters) + 1, column=0, columnspan=4, pady=10)

        # Create buttons
        cancel_button = tk.Button(
            button_frame,
            text="Cancel",
            command=lambda: self.on_button_click("Cancel"),
            width=20,
        )
        ok_button = tk.Button(
            button_frame,
            text="OK",
            command=lambda: self.on_button_click("OK"),
            width=20,
        )

        cancel_button.grid(row=0, column=0, padx=10)
        ok_button.grid(row=0, column=1, padx=10)

    def on_button_click(self, button_text):
        if button_text in ["Batch", "OK"]:
            self.win.destroy()
            # Handle Batch or OK button click
            print(f"{button_text} button clicked")

    def update(self):
        values = {}
        for param, (min_entry, max_entry, invert_var) in zip(
            self.parameters, self.entries
        ):
            values[param] = {
                "min": min_entry.get(),
                "max": max_entry.get(),
                "invert": invert_var.get(),
            }
        self.update_func(values)
        self.win.destroy()


class SelectPointsWindow:
    def __init__(self, parameters, update_func, plot_settings=None):
        self.win = tk.Toplevel()
        self.win.title("Selection criteria")
        self.update_func = update_func
        self.parameters = parameters
        self.plot_settings = plot_settings

        # Create header labels
        tk.Label(self.win, text="Channel", width=20, anchor="w").grid(
            row=0, column=0, padx=5, sticky="we"
        )
        tk.Label(self.win, text="Min", width=13, anchor="w").grid(
            row=0, column=1, padx=5, sticky="w"
        )
        tk.Label(self.win, text="Max", width=14, anchor="w").grid(
            row=0, column=2, padx=5, sticky="w"
        )
        tk.Label(self.win, text="Derivative", anchor="w", width=8).grid(
            row=0, column=3, padx=0, sticky="w"
        )
        tk.Label(self.win, text="Invert", width=8, anchor="w").grid(
            row=0, column=4, padx=0, sticky="w"
        )

        self.entries = []
        for i, param in enumerate(parameters):
            tk.Label(self.win, text=param, width=20).grid(
                row=i + 1, column=0, padx=5, pady=4, sticky="we"
            )

            min_entry = tk.Entry(self.win, width=13)
            min_entry.grid(row=i + 1, column=1, padx=5, pady=4, sticky="w")

            max_entry = tk.Entry(self.win, width=14)
            max_entry.grid(row=i + 1, column=2, padx=5, pady=4, sticky="w")

            derivative_var = tk.BooleanVar()
            derivative_check = tk.Checkbutton(self.win, variable=derivative_var)
            derivative_check.grid(row=i + 1, column=3, padx=8, pady=2, sticky="w")

            invert_var = tk.BooleanVar()
            invert_check = tk.Checkbutton(self.win, variable=invert_var)
            invert_check.grid(row=i + 1, column=4, padx=8, pady=2, sticky="w")

            self.entries.append((min_entry, max_entry, derivative_var, invert_var))

        # Create a custom progress bar using Canvas
        self.progress_canvas = tk.Canvas(
            self.win, width=550, height=10, bg="white", bd=0, highlightthickness=0
        )
        self.progress_canvas.grid(
            row=len(parameters) + 1, column=0, columnspan=5, pady=1
        )
        self.progress_bar = self.progress_canvas.create_rectangle(
            0, 0, 0, 10, fill="blue"
        )

        # Create a frame for the buttons
        button_frame = tk.Frame(self.win)
        button_frame.grid(row=len(parameters) + 2, column=0, columnspan=4, pady=10)

        # Create buttons
        cancel_button = tk.Button(
            button_frame,
            text="Cancel",
            command=lambda: self.on_button_click("Cancel"),
            width=20,
        )
        batch_button = tk.Button(
            button_frame,
            text="Batch",
            command=lambda: self.on_button_click("Batch"),
            width=20,
        )
        ok_button = tk.Button(
            button_frame,
            text="OK",
            command=lambda: self.on_button_click("OK"),
            width=20,
        )

        cancel_button.grid(row=0, column=0, padx=10)
        batch_button.grid(row=0, column=1, padx=10)
        ok_button.grid(row=0, column=2, padx=10)

    def on_button_click(self, button_text):
        def process_to_do(self, to_do):
            for i, (f, label) in enumerate(to_do):
                self.update_progress(100 - (i + 1) / len(to_do) * 100)
                self.update_func(f, label)

        if button_text in ["Cancel", "OK"]:
            self.win.destroy()
            # Handle Batch or OK button click
            print(f"{button_text} button clicked")
        elif button_text == "Batch":
            to_do = []
            for f in self.plot_settings["filenames"]:
                data = tio.hdf_data(f, "0")
                for label in data.list_labels():
                    to_do.append((f, label))
                    print(f, label)

    def update(self):
        values = {}
        for param, (min_entry, max_entry, invert_var) in zip(
            self.parameters, self.entries
        ):
            values[param] = {
                "min": min_entry.get(),
                "max": max_entry.get(),
                "invert": invert_var.get(),
            }
        self.update_func(values)
        self.win.destroy()


class SetFitParametersWindow:
    def __init__(self, parameters, update_func):
        self.win = tk.Toplevel()
        self.win.title("Set fit constraints")
        self.update_func = update_func
        self.parameters = parameters

        # Create header labels
        tk.Label(self.win, text="Parameter", width=20, anchor="w").grid(
            row=0, column=0, padx=5, sticky="we"
        )
        tk.Label(self.win, text="Value", width=13, anchor="w").grid(
            row=0, column=1, padx=5, sticky="w"
        )
        tk.Label(self.win, text="Min", width=13, anchor="w").grid(
            row=0, column=2, padx=5, sticky="w"
        )
        tk.Label(self.win, text="Max", width=14, anchor="w").grid(
            row=0, column=3, padx=5, sticky="w"
        )
        tk.Label(self.win, text="Type", width=5, anchor="w").grid(
            row=0, column=4, padx=0, sticky="w"
        )

        self.entries = []
        for i, param in enumerate(parameters):
            tk.Label(self.win, text=param, width=20).grid(
                row=i + 1, column=0, padx=5, pady=4, sticky="we"
            )

            value_entry = tk.Entry(self.win, width=13)
            value_entry.grid(row=i + 1, column=1, padx=5, pady=4, sticky="w")

            min_entry = tk.Entry(self.win, width=13)
            min_entry.grid(row=i + 1, column=2, padx=5, pady=4, sticky="w")

            max_entry = tk.Entry(self.win, width=14)
            max_entry.grid(row=i + 1, column=3, padx=5, pady=4, sticky="w")

            type_var = tk.StringVar()
            type_box = ttk.Combobox(
                self.win,
                values=["Fit", "Global", "Local"],
                width=10,
                textvariable=type_var,
            )
            type_box.grid(row=i + 1, column=4, padx=5, pady=4, sticky="w")

            self.entries.append((value_entry, min_entry, max_entry, type_var))

        # Create a frame for the buttons
        button_frame = tk.Frame(self.win)
        button_frame.grid(row=len(parameters) + 1, column=0, columnspan=4, pady=10)

        # Create buttons
        cancel_button = tk.Button(
            button_frame,
            text="Cancel",
            command=lambda: self.on_button_click("Cancel"),
            width=20,
        )
        batch_button = tk.Button(
            button_frame,
            text="Batch",
            command=lambda: self.on_button_click("Batch"),
            width=20,
        )
        ok_button = tk.Button(
            button_frame,
            text="OK",
            command=lambda: self.on_button_click("OK"),
            width=20,
        )

        cancel_button.grid(row=0, column=0, padx=10)
        batch_button.grid(row=0, column=1, padx=10)
        ok_button.grid(row=0, column=2, padx=10)

    def on_button_click(self, button_text):
        if button_text in ["Cancel", "Batch", "OK"]:
            self.win.destroy()
            # Handle Batch or OK button click
            print(f"{button_text} button clicked")

    def update(self):
        values = {}
        for param, (min_entry, max_entry, invert_var) in zip(
            self.parameters, self.entries
        ):
            values[param] = {
                "min": min_entry.get(),
                "max": max_entry.get(),
                "invert": invert_var.get(),
            }
        self.update_func(values)
        self.win.destroy()


# class FitParametersWindow:
#     def __init__(self, parameters, update_func):
#         self.win = tk.Toplevel()
#         self.win.title("Selection criteria")
#         self.update_func = update_func

#         # Create a frame for the header labels
#         header_frame = tk.Frame(self.win)
#         header_frame.grid(row=0, column=0, columnspan=4, pady=10)

#         # Create header labels
#         tk.Label(header_frame, text="Channel", width=10).grid(
#             row=0, column=0, padx=5, sticky="w")
#         tk.Label(header_frame, text="Min", width=12).grid(
#             row=0, column=1, padx=10, sticky="w")
#         tk.Label(header_frame, text="Max", width=12).grid(
#             row=0, column=2, padx=10, sticky="w")
#         tk.Label(header_frame, text="     Invert", width=10).grid(
#             row=0, column=3, padx=1, sticky="w")
#         tk.Label(header_frame, text="Derivative", width=10).grid(
#             row=0, column=4, padx=1, sticky="w")


#         # Create a frame for the parameter rows
#         param_frame = tk.Frame(self.win)
#         param_frame.grid(row=1, column=0, columnspan=4)

#         # Create parameter rows
#         entry_vars = {}
#         for i, param in enumerate(parameters):
#             tk.Label(param_frame, text=param, width=15).grid(
#                 row=i, column=0, padx=5, pady=0, sticky="e")

#             value_var = tk.StringVar()
#             min_var = tk.StringVar()
#             max_var = tk.StringVar()
#             invert_var = tk.BooleanVar()
#             derivative_var = tk.BooleanVar

#             value_entry = tk.Entry(param_frame, textvariable=value_var, width=15)
#             min_entry = tk.Entry(param_frame, textvariable=min_var, width=15)
#             max_entry = tk.Entry(param_frame, textvariable=max_var, width=15)
#             invert_check = tk.Checkbutton(param_frame, variable=invert_var)
#             derivative_check = tk.Checkbutton(param_frame, variable=derivative_var)

#             value_entry.grid(row=i, column=1, padx=5, pady=0, sticky="w")
#             min_entry.grid(row=i, column=2, padx=5, pady=0, sticky="w")
#             max_entry.grid(row=i, column=3, padx=5, pady=0, sticky="w")
#             invert_check.grid(row=i, column=4, padx=15, pady=0, sticky="w")
#             derivative_check.grid(row=i, column=5, padx=15, pady=0, sticky="w")

#             entry_vars[param] = {"value": value_var,"min": min_var,
#                                 "max": max_var, "invert": invert_var}

#             # Create a frame for the buttons
#             button_frame = tk.Frame(self.win)
#             button_frame.grid(row=2, column=0, columnspan=4, pady=10)

#             # Create buttons
#             cancel_button = tk.Button(button_frame, text="Cancel",
#                                     command=lambda: on_button_click("Cancel"), width=15)
#             batch_button = tk.Button(button_frame, text="Batch",
#                                 command=lambda: on_button_click("Batch"), width=15)
#             ok_button = tk.Button(button_frame, text="OK",
#                                 command=lambda: on_button_click("OK"), width=15)

#             cancel_button.grid(row=0, column=0, padx=10)
#             batch_button.grid(row=0, column=2, padx=10)
#             ok_button.grid(row=0, column=3, padx=10)

#         def on_button_click(button_text):
#             if button_text in ["Batch", "OK"]:
#                 self.win.destroy()
#                 # Handle Batch or OK button click
#                 print(f"{button_text} button clicked")
#                 return_values = {param: {"min": vars["min"].get(), "max": vars["max"].get(
#                 ), "invert": vars["invert"].get()} for param, vars in entry_vars.items()
#                 if vars["min"].get() != "" or vars["max"].get() != ""}
#                 self.update_func(return_values)
#             elif button_text == "Cancel":
#                 self.win.destroy()
#                 # Handle Cancel button click
#                 print("Cancel button clicked")


class MainWindow:
    def __init__(self, master):
        self.master = master
        self.settings = None
        self.SettingsValue = None

        self.frame = tk.Frame(self.master, width=200, height=200)
        self.frame.pack()

        self.button = tk.Button(
            self.frame, text="Settings Window", command=self.openSettings
        )
        self.button.place(x=50, y=50)

        self.button2 = tk.Button(
            self.frame, text="Update Settings", command=self.updateLabelinSettings
        )
        self.button2.place(x=50, y=150)

    def update(self, n):
        self.SettingsValue = n
        print(self.SettingsValue)

    def openSettings(self):
        parameters = ["Parameter1", "Parameter2", "Parameter3", "Parameter4"]
        self.settings = SelectPointsWindow(parameters, self.update)

    def updateLabelinSettings(self):
        prompt = askinteger("Input", "Input an Integer")
        if self.settings != None:
            self.settings.label.configure(text="Updated")


if __name__ == "__main__":
    win = tk.Tk()
    window = MainWindow(win)
    win.mainloop()
