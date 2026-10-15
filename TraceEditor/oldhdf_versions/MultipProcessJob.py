import time
import multiprocessing
import keyboard
import numpy as np
from icecream import ic
import pandas as pd

import TraceIO as tio
import ForceSpectroscopy as fs


class CalculationProcess(multiprocessing.Process):
    def __init__(self, jobs, parameters=None, fitmodels=None):
        super(CalculationProcess, self).__init__()
        self.result_queue = multiprocessing.Queue()
        self.stop_event = multiprocessing.Event()
        self.jobs = jobs
        self.parameters = parameters
        self.fitmodels = fitmodels
        self.start()

    def run(self):
        for i, job in enumerate(self.jobs):
            if self.stop_event.is_set():
                break
            data = tio.hdf_data(job["filename"], job["label"])
            # Do calculation
            if job["task"] == "select points":
                data.select_points(self.parameters)
            elif job["task"] == "select traces":
                data.select_trace(self.parameters)
            elif job["task"] in self.fitmodels if self.fitmodels else []:
                data.fit(job["task"], self.parameters)
            elif job["task"] == "Calc free DNA": 
                data.calc_free_DNA()
            elif job["task"] == "calibrate force":
                data.calibrate_force(job["parameters"])
            elif job["task"] == "zMagnetForce":
                data.recalculate_force(job["parameters"])
            else:
                time.sleep(0.1)
            self.result_queue.put([job])
        self.result_queue.put([])

    def stop(self):
        self.result_queue.put([])
        self.stop_event.set()
        self.join()
        print("Process stopped")


def main():
    jobs = [
        ("none", "C:/tmp/vsc/TraceEditor/data/data_014.hdf", "0"),
        ("none", "C:/tmp/vsc/TraceEditor/data/data_014.hdf", "1"),
    ]
    process = CalculationProcess(jobs)
    process.start()
    keyboard.add_hotkey("esc", process.stop)

    while True:
        message = process.result_queue.get()
        if message == "DONE":
            break
        else:
            print(message)

    process.stop()


if __name__ == "__main__":
    main()
