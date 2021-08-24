import sys, yaml, json, random, os, re
from matplotlib.backends.backend_qt5 import NavigationToolbar2QT
import threading, faulthandler
import gc as garbage

faulthandler.enable()

# Import QApplication and the required widgets from PyQt5.QtWidgets
from PyQt5.QtCore import *
from PyQt5.QtWidgets import *
from PyQt5.QtGui import *

# Import matplotlib for heatmap rendering
import matplotlib
import networkx as nx

matplotlib.use("QT5Agg")
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
import matplotlib.pyplot as plt

# Import pandas and numpy for data manipulation
import numpy as np
import pandas as pd

# Import mfold required libraries
from mfold_library import Region, Mfold
from genetic import Sequence

# abl20 able25, abre25 abr20, ABRE25 abd36 ABLE25, able25 ABD36 abre25, acl20 acle25, acre25 acr20, ACRE25 acd36 ACLE25, acle25 ACD36 acre25
# able25 ABB36 abre25, ABRE25 abb36 ABLE25, acle25 ACB36 acre25, ACRE25 acb36 ACLE25, adle25 ADB36 adre25, ADRE25 adb36 ADLE25

# To do:
# output the 2D shape of the input given V
# Label plot numbers and strands V
# fix the progress bar update in the loop calculation V
# Generate figures only for the highlighted regions V
# tkinter after: V
def clear_files():
    extensions = [
        ".seq",
        ".ct",
        ".det",
        ".out",
        ".pnt",
        ".sav",
        ".test",
        ".log",
        ".cmd",
        ".aux",
        ".con",
    ]
    cwd = os.getcwd()
    test = os.listdir(cwd)

    for item in test:
        for extension in extensions:
            if item.endswith(extension):
                os.remove(os.path.join(cwd, item))


class Randomizer(QObject):
    render = pyqtSignal()
    finished = pyqtSignal()
    error = pyqtSignal()
    progress = pyqtSignal(str)
    highlight = pyqtSignal()
    progression = pyqtSignal(int)
    trials = 0
    randomizations = 0

    def __init__(self, parent):
        QObject.__init__(self)
        self.parent = parent

    @pyqtSlot()
    def run_randomization(self):
        params = {}
        params["raw_structure"] = self.parent.raw_structure
        params["energy_matrix"] = self.parent.energy
        params["fixed_regions"] = self.parent.fixed_regions
        params["input_sequence_definitions"] = self.parent.input_sequence
        params["highlighted"] = self.parent.highlighted
        params["index"] = self.parent.index
        params["shape"] = self.parent.shape

        self.progression.emit(int((self.trials / 5) * 100))
        new = 0
        max = 0
        for iteration in range(6):
            flag = False
            counter = 0
            self.progress.emit("Randomizing")
            while counter < len(self.parent.highlighted):
                fields = []

                for a in range(len(self.parent.highlighted)):
                    text = ""

                    for b in range(len(self.parent.highlighted[a])):
                        if self.parent.highlighted[a][b] == "]":
                            flag = False

                        elif flag:
                            text += random.choice(["G", "T", "C", "A"])

                        elif self.parent.highlighted[a][b] == "[":
                            flag = True

                        else:
                            text += self.parent.highlighted[a][b]

                    fields.append(text)

                gc = []
                while len(self.parent.population_size) != len(self.parent.header):
                    _ = 0
                data = gc_content(gc, self)
                gci = 0
                for _ in fields[counter]:
                    if _ == "C" or _ == "G":
                        gci += 1

                if gc[self.parent.index[counter]] != gci:
                    if "]" in self.parent.highlighted[counter]:
                        flag = False

                    if (
                        "[" in self.parent.highlighted[counter]
                        and not "]" in self.parent.highlighted[counter]
                    ):
                        flag = True
                    self.randomizations += 1
                    if self.randomizations > 20:
                        while self.parent.center_layout.count() != 0:
                            self.parent.center_layout.removeRow(0)

                        self.parent.header.clear()
                        self.parent.field.clear()
                        self.parent.label.clear()
                        self.parent.check.clear()

                        self.parent.raw_structure = params["raw_structure"]
                        self.parent.energy = params["energy_matrix"]
                        self.parent.fixed_regions = params["fixed_regions"]
                        self.parent.input_sequence = params[
                            "input_sequence_definitions"
                        ]
                        self.parent.index = params["index"]
                        self.parent.highlighted = params["highlighted"]
                        self.parent.shape = params["shape"]

                        self.render.emit()
                        self.randomizations = 0
                        garbage.collect()
                        break
                    else:
                        continue

                else:
                    self.randomizations = 0
                    if (
                        self.parent.header[self.parent.index[counter]]
                        in self.parent.fixed_regions.keys()
                    ):
                        error = QErrorMessage(self)
                        error.showMessage(
                            f"Please make sure to uncheck the field {self.parent.header[self.parent.index[counter]].lower()} to allow modifications"
                        )
                        break

                    if self.parent.index.count(self.parent.index[counter]) > 1:
                        indices = [
                            i
                            for i, x in enumerate(self.parent.index)
                            if x == self.parent.index[counter]
                        ]
                        indices.pop(0)

                        for i in indices:
                            temp = self.parent.highlighted[i]

                            new_indices_o = [
                                j for j in range(len(temp)) if temp.startswith("[", j)
                            ]

                            new_indices_c = [
                                j for j in range(len(temp)) if temp.startswith("]", j)
                            ]
                            temp = [c for c in fields[counter]]
                            for k in new_indices_o:
                                temp.insert(k, "[")

                            for l in new_indices_c:
                                temp.insert(l, "]")

                            temp = "".join(temp)
                            self.parent.highlighted[i] = temp

                    for i in range(len(self.parent.header)):
                        if self.parent.header[self.parent.index[counter]].islower():
                            if (
                                self.parent.header[i]
                                == self.parent.header[
                                    self.parent.index[counter]
                                ].upper()
                            ):
                                indices = [
                                    j for j, x in enumerate(self.parent.index) if x == i
                                ]
                                for k in indices:
                                    temp = self.parent.highlighted[k]
                                    new_indices_o = [
                                        j
                                        for j in range(len(temp))
                                        if temp.startswith("[", j)
                                    ]

                                    new_indices_c = [
                                        j
                                        for j in range(len(temp))
                                        if temp.startswith("]", j)
                                    ]
                                    temp = (
                                        fields[counter][::-1]
                                        .replace("A", "temp")
                                        .replace("T", "A")
                                        .replace("temp", "T")
                                    )
                                    temp = (
                                        temp.replace("C", "temp")
                                        .replace("G", "C")
                                        .replace("temp", "G")
                                    )
                                    temp = [c for c in temp]
                                    for l in new_indices_o:
                                        temp.insert(l, "[")

                                    for m in new_indices_c:
                                        temp.insert(m, "]")

                                    temp = "".join(temp)
                                    self.parent.highlighted[k] = temp

                        if self.parent.header[self.parent.index[counter]].isupper():
                            if (
                                self.parent.header[i]
                                == self.parent.header[
                                    self.parent.index[counter]
                                ].lower()
                            ):
                                indices = [
                                    j for j, x in enumerate(self.parent.index) if x == i
                                ]
                                for k in indices:
                                    temp = self.parent.highlighted[k]

                                    new_indices_o = [
                                        j
                                        for j in range(len(temp))
                                        if temp.startswith("[", j)
                                    ]

                                    new_indices_c = [
                                        j
                                        for j in range(len(temp))
                                        if temp.startswith("]", j)
                                    ]
                                    temp = (
                                        fields[counter][::-1]
                                        .replace("A", "temp")
                                        .replace("T", "A")
                                        .replace("temp", "T")
                                    )
                                    temp = (
                                        temp.replace("C", "temp")
                                        .replace("G", "C")
                                        .replace("temp", "G")
                                    )
                                    temp = [c for c in temp]
                                    for l in new_indices_o:
                                        temp.insert(l, "[")
                                    for m in new_indices_c:
                                        temp.insert(m, "]")

                                    temp = "".join(temp)

                                    self.parent.highlighted[k] = temp

                    self.parent.field[self.parent.index[counter]].setText(
                        fields[counter]
                    )
                    self.parent.field[self.parent.index[counter]].setModified(True)
                    self.parent.strand_update()

                    try:
                        if "]" in self.parent.highlighted[counter]:
                            if "[" in self.parent.highlighted[counter]:
                                counter += 1
                                continue
                            flag = True
                    except IndexError:
                        pass

                    counter += 1

            higher = 0
            for i in range(len(self.parent.energy)):
                for j in range(len(self.parent.energy)):
                    if self.parent.energy[i][j] != None:
                        if self.parent.energy[i][j] > higher:
                            higher = self.parent.energy[i][j]

            max = higher

            self.progress.emit("Testing randomization")

            self.parent.btn3.setEnabled(False)
            mfold = Mfold(output_folder="./", mfold_command="mfold_quik")

            self.parent.energy = [
                [None for strand1 in self.parent.strand]
                for strand2 in self.parent.strand
            ]

            for i, strand1 in enumerate(self.parent.strand):
                for j, strand2 in enumerate(self.parent.strand):
                    mfold.run(strand1, strand2, f"{i}_{j}.seq", f"{i}_{j}.aux")
                    try:
                        with open(f"{i}_{j}.det", "r") as configfile:
                            for line in configfile:
                                if line.startswith(" dG = "):
                                    self.parent.energy[i][j] = abs(float(line[10:15]))
                                    break
                    except FileNotFoundError:
                        self.parent.energy[i][j] = 0

            higher = 0
            for i in range(len(self.parent.energy)):
                for j in range(len(self.parent.energy)):
                    if self.parent.energy[i][j] != None:
                        if self.parent.energy[i][j] > higher:
                            higher = self.parent.energy[i][j]

            new = higher
            if new > max and self.trials < 5:
                while self.parent.center_layout.count() != 0:
                    self.parent.center_layout.removeRow(0)

                self.parent.header.clear()
                self.parent.field.clear()
                self.parent.label.clear()
                self.parent.check.clear()

                self.parent.raw_structure = params["raw_structure"]
                self.parent.energy = params["energy_matrix"]
                self.parent.fixed_regions = params["fixed_regions"]
                self.parent.input_sequence = params["input_sequence_definitions"]
                self.parent.highlighted = params["highlighted"]
                self.parent.index = params["index"]
                self.parent.shape = params["shape"]

                self.render.emit()
                clear_files()

                self.trials += 1
                self.progression.emit(int((self.trials / 5) * 100))

            else:
                break

        flag = False
        if new > max:
            while self.parent.center_layout.count() != 0:
                self.parent.center_layout.removeRow(0)

            self.parent.header.clear()
            self.parent.field.clear()
            self.parent.label.clear()
            self.parent.check.clear()

            self.parent.raw_structure = params["raw_structure"]
            self.parent.energy = params["energy_matrix"]
            self.parent.fixed_regions = params["fixed_regions"]
            self.parent.input_sequence = params["input_sequence_definitions"]
            self.parent.index = params["index"]
            self.parent.highlighted = params["highlighted"]
            self.parent.shape = params["shape"]

            self.render.emit()
            flag = True

        self.progress.emit("Generating structures")
        self.progression.emit(100)
        self.parent.automate = False

        find_thief(self)
        self.highlight.emit()

        self.progress.emit("DONE")

        if flag:
            if not self.parent.automate:
                self.error.emit()

        self.parent.btn3.setEnabled(True)
        self.finished.emit()


class Looper(QObject):
    """
    A class that is used as a worker/thread for loop calculations. It is defined by signals that are sent to the main thread and updates
    the main self class

    Args:
            (from main class)
            self.btn3 : "calculate" button
            self.automate : flag responsable for loop calculation
            self.file_counter : counter of files created and the current loop calculation
            self.iterator : counter for max iterations of loop calculation
    """

    finished = pyqtSignal()
    highlight = pyqtSignal()
    progression = pyqtSignal(str)
    render = pyqtSignal()
    progress = pyqtSignal(int)
    randomizations = 0
    trials = 0

    def __init__(self, parent):
        QObject.__init__(self)
        self.parent = parent

    @pyqtSlot()
    def run_routine(self):
        """ "
        This function calls the mfold thread and waits for it to finish before updating the progress bar and sends the progress to
        the main class.

        Args:
                self.progression : responsable for
                self.progress :
                self.finished :
        """
        self.parent.btn3.setEnabled(False)
        self.parent.btn8.setEnabled(False)
        for count in range(self.parent.iterator):
            self.progress.emit(
                int((self.parent.file_counter / self.parent.iterator) * 100)
            )
            clear_files()

            self.parent.automate = True
            self.progression.emit("Calculating")
            mfold = Mfold(output_folder="./", mfold_command="mfold_quik")

            self.parent.energy = [
                [None for strand1 in self.parent.strand]
                for strand2 in self.parent.strand
            ]

            for a, strand1 in enumerate(self.parent.strand):
                for b, strand2 in enumerate(self.parent.strand):
                    mfold.run(strand1, strand2, f"{a}_{b}.seq", f"{a}_{b}.aux")
                    try:
                        with open(f"{a}_{b}.det", "r") as configfile:
                            for line in configfile:
                                if line.startswith(" dG = "):
                                    self.parent.energy[a][b] = abs(float(line[10:15]))
                                    break
                    except FileNotFoundError:
                        self.parent.energy[a][b] = 0

            # garbage.collect()
            find_thief(self)
            self.parent.highlighted.clear()
            self.parent.index.clear()
            self.parent.shape.clear()

            for g in range(2):
                intermediate_process(self.parent)

            self.progression.emit("Randomizing")

            params = {}
            data = gc_content([], self)
            for _ in range(len(self.parent.header)):
                self.parent.input_sequence[self.parent.header[_]] = data[
                    self.parent.header[_]
                ]

            del data
            # garbage.collect()
            params["raw_structure"] = self.parent.raw_structure
            params["fixed_regions"] = self.parent.fixed_regions
            params["input_sequence_definitions"] = self.parent.input_sequence
            params["energy_matrix"] = self.parent.energy
            params["highlighted"] = self.parent.highlighted
            params["index"] = self.parent.index
            params["shape"] = self.parent.shape

            new = 0
            max = 0

            for iteration in range(6):
                flag = False
                counter = 0
                while counter < len(self.parent.highlighted):
                    fields = []

                    for a in range(len(self.parent.highlighted)):
                        text = ""

                        for b in range(len(self.parent.highlighted[a])):
                            if self.parent.highlighted[a][b] == "]":
                                flag = False

                            elif flag:
                                text += random.choice(["G", "T", "C", "A"])

                            elif self.parent.highlighted[a][b] == "[":
                                flag = True

                            else:
                                text += self.parent.highlighted[a][b]

                        fields.append(text)

                    gc = []
                    while len(self.parent.population_size) != len(self.parent.header):
                        _ = 0
                    temporary = gc_content(gc, self)
                    gci = 0
                    for _ in fields[counter]:
                        if _ == "C" or _ == "G":
                            gci += 1

                    if gc[self.parent.index[counter]] != gci:
                        if "]" in self.parent.highlighted[counter]:
                            flag = False

                        if (
                            "[" in self.parent.highlighted[counter]
                            and not "]" in self.parent.highlighted[counter]
                        ):
                            flag = True
                        self.randomizations += 1
                        if self.randomizations > 20:
                            while self.parent.center_layout.count() != 0:
                                self.parent.center_layout.removeRow(0)

                            self.parent.header.clear()
                            self.parent.field.clear()
                            self.parent.label.clear()
                            self.parent.check.clear()
                            self.parent.input_sequence = {}

                            self.parent.raw_structure = params["raw_structure"]
                            self.parent.energy = params["energy_matrix"]
                            self.parent.fixed_regions = params["fixed_regions"]
                            self.parent.input_sequence = params[
                                "input_sequence_definitions"
                            ]
                            self.parent.highlighted = params["highlighted"]
                            self.parent.index = params["index"]
                            self.parent.shape = params["shape"]

                            self.render.emit()

                            self.randomizations = 0
                            # garbage.collect()
                            break
                        else:
                            continue

                    else:
                        # garbage.collect()
                        self.randomizations = 0
                        if (
                            self.parent.header[self.parent.index[counter]]
                            in self.parent.fixed_regions.keys()
                        ):
                            error = QErrorMessage(self)
                            error.showMessage(
                                f"Please make sure to uncheck the field {self.parent.header[self.parent.index[counter]].lower()} to allow modifications"
                            )
                            break

                        if self.parent.index.count(self.parent.index[counter]) > 1:
                            indices = [
                                i
                                for i, x in enumerate(self.parent.index)
                                if x == self.parent.index[counter]
                            ]
                            indices.pop(0)

                            for i in indices:
                                temp = self.parent.highlighted[i]

                                new_indices_o = [
                                    j
                                    for j in range(len(temp))
                                    if temp.startswith("[", j)
                                ]

                                new_indices_c = [
                                    j
                                    for j in range(len(temp))
                                    if temp.startswith("]", j)
                                ]
                                temp = [c for c in fields[counter]]
                                for k in new_indices_o:
                                    temp.insert(k, "[")

                                for l in new_indices_c:
                                    temp.insert(l, "]")

                                temp = "".join(temp)
                                self.parent.highlighted[i] = temp

                        for i in range(len(self.parent.header)):
                            if self.parent.header[self.parent.index[counter]].islower():
                                if (
                                    self.parent.header[i]
                                    == self.parent.header[
                                        self.parent.index[counter]
                                    ].upper()
                                ):
                                    indices = [
                                        j
                                        for j, x in enumerate(self.parent.index)
                                        if x == i
                                    ]
                                    for _ in indices:
                                        temp = self.parent.highlighted[_]
                                        new_indices_o = [
                                            j
                                            for j in range(len(temp))
                                            if temp.startswith("[", j)
                                        ]

                                        new_indices_c = [
                                            j
                                            for j in range(len(temp))
                                            if temp.startswith("]", j)
                                        ]
                                        temp = (
                                            fields[counter][::-1]
                                            .replace("A", "temp")
                                            .replace("T", "A")
                                            .replace("temp", "T")
                                        )
                                        temp = (
                                            temp.replace("C", "temp")
                                            .replace("G", "C")
                                            .replace("temp", "G")
                                        )
                                        temp = [c for c in temp]
                                        for l in new_indices_o:
                                            temp.insert(l, "[")

                                        for m in new_indices_c:
                                            temp.insert(m, "]")

                                        temp = "".join(temp)
                                        self.parent.highlighted[_] = temp

                            if self.parent.header[self.parent.index[counter]].isupper():
                                if (
                                    self.parent.header[i]
                                    == self.parent.header[
                                        self.parent.index[counter]
                                    ].lower()
                                ):
                                    indices = [
                                        j
                                        for j, x in enumerate(self.parent.index)
                                        if x == i
                                    ]
                                    for k in indices:
                                        temp = self.parent.highlighted[k]

                                        new_indices_o = [
                                            j
                                            for j in range(len(temp))
                                            if temp.startswith("[", j)
                                        ]

                                        new_indices_c = [
                                            j
                                            for j in range(len(temp))
                                            if temp.startswith("]", j)
                                        ]
                                        temp = (
                                            fields[counter][::-1]
                                            .replace("A", "temp")
                                            .replace("T", "A")
                                            .replace("temp", "T")
                                        )
                                        temp = (
                                            temp.replace("C", "temp")
                                            .replace("G", "C")
                                            .replace("temp", "G")
                                        )
                                        temp = [c for c in temp]
                                        for l in new_indices_o:
                                            temp.insert(l, "[")
                                        for m in new_indices_c:
                                            temp.insert(m, "]")

                                        temp = "".join(temp)

                                        self.parent.highlighted[k] = temp

                        self.parent.field[self.parent.index[counter]].setText(
                            fields[counter]
                        )
                        self.parent.field[self.parent.index[counter]].setModified(True)
                        self.parent.strand_update()

                        try:
                            if "]" in self.parent.highlighted[counter]:
                                if "[" in self.parent.highlighted[counter]:
                                    counter += 1
                                    continue
                                flag = True
                        except IndexError:
                            pass

                        counter += 1

                higher = 0
                for i in range(len(self.parent.energy)):
                    for j in range(len(self.parent.energy)):
                        if self.parent.energy[i][j] != None:
                            if self.parent.energy[i][j] > higher:
                                higher = self.parent.energy[i][j]

                max = higher

                mfold = Mfold(output_folder="./", mfold_command="mfold_quik")

                self.parent.energy = [
                    [None for strand1 in self.parent.strand]
                    for strand2 in self.parent.strand
                ]

                for i, strand1 in enumerate(self.parent.strand):
                    for j, strand2 in enumerate(self.parent.strand):
                        mfold.run(strand1, strand2, f"{i}_{j}.seq", f"{i}_{j}.aux")
                        try:
                            with open(f"{i}_{j}.det", "r") as configfile:
                                for line in configfile:
                                    if line.startswith(" dG = "):
                                        self.parent.energy[i][j] = abs(
                                            float(line[10:15])
                                        )
                                        break
                        except FileNotFoundError:
                            self.parent.energy[i][j] = 0

                higher = 0
                for i in range(len(self.parent.energy)):
                    for j in range(len(self.parent.energy)):
                        if self.parent.energy[i][j] != None:
                            if self.parent.energy[i][j] > higher:
                                higher = self.parent.energy[i][j]

                new = higher
                if new >= max and self.trials < 5:
                    while self.parent.center_layout.count() != 0:
                        self.parent.center_layout.removeRow(0)

                    self.parent.header.clear()
                    self.parent.field.clear()
                    self.parent.label.clear()
                    self.parent.check.clear()
                    self.parent.input_sequence = {}

                    self.parent.raw_structure = params["raw_structure"]
                    self.parent.energy = params["energy_matrix"]
                    self.parent.fixed_regions = params["fixed_regions"]
                    self.parent.input_sequence = params["input_sequence_definitions"]
                    self.parent.highlighted = params["highlighted"]
                    self.parent.index = params["index"]
                    self.parent.shape = params["shape"]

                    self.render.emit()

                    self.trials += 1

                else:
                    break

            if new >= max:
                while self.parent.center_layout.count() != 0:
                    self.parent.center_layout.removeRow(0)

                self.parent.header.clear()
                self.parent.field.clear()
                self.parent.label.clear()
                self.parent.check.clear()
                self.parent.input_sequence = {}
                self.parent.raw_structure = params["raw_structure"]
                self.parent.energy = params["energy_matrix"]
                self.parent.fixed_regions = params["fixed_regions"]
                self.parent.input_sequence = params["input_sequence_definitions"]
                self.parent.highlighted = params["highlighted"]
                self.parent.index = params["index"]
                self.parent.shape = params["shape"]
                self.render.emit()

                flag = True

            # find_thief(self)
            # intermediate_process(self.parent)

            # garbage.collect()

            self.parent.file_counter += 1
            garbage.collect()
            # print("here")
            # print(new, max)

        find_thief(self)
        self.parent.highlighted.clear()
        self.parent.index.clear()
        self.parent.shape.clear()
        self.highlight.emit()
        self.progress.emit(100)
        self.parent.automate = False
        self.progression.emit("DONE")
        self.parent.file_counter = 0
        self.parent.iterator = 0
        self.finished.emit()
        self.parent.btn8.setEnabled(True)
        self.parent.btn3.setEnabled(True)


class Worker(QObject):
    """
    A class that is used as a worker/thread to call the mfold program and update the energy matrix in the main class.

    Args:
            mfold : call to Mfold class to interface with the Mfold program
            (from main class)
            self.energy : store the energy loss returned by the Mfold program
            self.automate : flag responsable to control the progress bar depending on the loop calculation or a single calculation
    """

    finished = pyqtSignal()
    progress = pyqtSignal(str)
    render = pyqtSignal()
    progression = pyqtSignal(int)

    def __init__(self, parent):
        QObject.__init__(self)
        self.parent = parent

    @pyqtSlot()
    def run_call(self):
        """
        This function initiates the call to the mfold program and populates the energy matrix of the main class while updating the progress.

        Args:
                self.finished :
                self.progress :
                self.progression :
        """
        mfold = Mfold(output_folder="./", mfold_command="mfold_quik")

        self.progress.emit("PROGRESSING")

        self.parent.energy = [
            [None for strand1 in self.parent.strand] for strand2 in self.parent.strand
        ]

        for i, strand1 in enumerate(self.parent.strand):
            if not self.parent.automate:
                self.progression.emit(int((i / len(self.parent.strand)) * 100))

            for j, strand2 in enumerate(self.parent.strand):
                mfold.run(strand1, strand2, f"{i}_{j}.seq", f"{i}_{j}.aux")
                try:
                    with open(f"{i}_{j}.det", "r") as configfile:
                        for line in configfile:
                            if line.startswith(" dG = "):
                                self.parent.energy[i][j] = abs(float(line[10:15]))
                                break
                except FileNotFoundError:
                    self.parent.energy[i][j] = 0

        if not self.parent.automate:
            self.progression.emit(100)

        find_thief(self)
        self.render.emit()
        self.finished.emit()

        self.progress.emit("DONE")
        garbage.collect()


def find_thief(self):
    """
    This function creates a dataframe by reading the ".ct" files generated by mfold and storing information relative to desired regions.

    Args:
            (from main class)
            self.strand : Structure class used to define the strand sequences
            self.regions : Region class used to define the regions defined
            self.energy : Matrix used to store the energy loss returned by the Mfold program
            self.thief : Dataframe containing the data generated by Mfold regarding the specific regions responsable for the energy loss
    """
    data = {}
    flag = False
    for i in range(len(self.parent.strand)):
        for j in range(len(self.parent.strand)):
            data[f"{i}{j}"] = {}

            # In case of already complementary regions
            list1 = self.parent.regions[i].lower().split("--")
            list2 = self.parent.regions[j].lower().split("--")
            list2 = list2[::-1]

            if list1 == list2 or self.parent.energy[i][j] == 0:
                self.parent.energy[i][j] = 0
                data[f"{i}{j}"]["energy"] = 0
                data[f"{i}{j}"]["letter"] = None
                data[f"{i}{j}"]["max"] = 0
                data[f"{i}{j}"]["index_app"] = 0
                data[f"{i}{j}"]["index_ct"] = 0
                data[f"{i}{j}"]["end"] = 0
                data[f"{i}{j}"]["shape"] = None

            else:
                try:
                    with open(f"{i}_{j}.ct", "r") as configfile:
                        counter = 0
                        max = 0
                        for line in configfile:
                            try:
                                if flag and int(line[26:30]) == 0:
                                    if counter <= len(self.parent.strand[i].bases):
                                        end = counter
                                    else:
                                        end = counter - 3
                                    flag = False
                                if flag:
                                    letter += line[7:8]

                                if int(line[26:30]) > max:
                                    letter = ""
                                    max = int(line[26:30])
                                    flag = True
                                    letter += line[7:8]

                                    if counter <= len(self.parent.strand[i].bases):
                                        index = counter
                                    else:
                                        index = counter - 3

                                    index1 = counter

                                if counter > (
                                    len(self.parent.strand[i].bases)
                                    + len(self.parent.strand[j].bases)
                                    + 4
                                ):
                                    break
                                counter += 1
                            except ValueError:
                                counter += 1
                                pass

                        data[f"{i}{j}"]["energy"] = self.parent.energy[i][j]
                        data[f"{i}{j}"]["letter"] = letter
                        data[f"{i}{j}"]["max"] = max
                        data[f"{i}{j}"]["index_app"] = index
                        data[f"{i}{j}"]["index_ct"] = index1
                        data[f"{i}{j}"]["end"] = end
                        data[f"{i}{j}"]["shape"] = f"{i}_{j}.png"

                except FileNotFoundError:
                    data[f"{i}{j}"]["energy"] = 0
                    data[f"{i}{j}"]["letter"] = None
                    data[f"{i}{j}"]["max"] = 0
                    data[f"{i}{j}"]["index_app"] = 0
                    data[f"{i}{j}"]["index_ct"] = 0
                    data[f"{i}{j}"]["end"] = 0
                    data[f"{i}{j}"]["shape"] = None
                    pass

    df = pd.DataFrame.from_dict(data)
    self.parent.thief = df


def intermediate_process(self):
    """
    This function is used to manipulate the dataframe generated by find_thief() and highlights the specific fields that are responsable for
    the energy loss.

    Args:
            (From main class)
            self.thief : Dataframe containing the data generated by Mfold regarding the specific regions responsable for the energy loss
            self.strand : Structure class used to define the strand sequences
            self.regions : Region class used to define the regions defined
            self.header : Region identifiers provided by the user in user_input()
            self.index : Indices of the fields of interest
            self.highlighted : Fields of interest to be highlighted
            self.field : Fields generated by the user in user_input()
    """
    if self.parent != None:
        self.parent = self
    temp = self.thief.transpose()

    for a in range(len(temp["energy"])):
        header = []
        if temp["energy"][a] == temp["energy"].max():
            end = temp["end"][a]
            atindex = temp["index_app"][a]

            _ = len(self.strand[int(self.thief.columns[a][0])].bases)
            if temp["end"][a] - 1 <= _:
                strand1 = int(self.thief.columns[a][0])
            else:
                strand1 = int(self.thief.columns[a][1])
                end = end - len(self.strand[int(self.thief.columns[a][0])].bases)
                atindex = atindex - len(
                    self.strand[int(self.thief.columns[a][0])].bases
                )

            for b in range(len(self.regions[strand1].split("--"))):
                header.append(self.regions[strand1].split("--")[b])

            data = gc_content([], self)
            length = [len(data[header[_]]) for _ in range(len(header))]

            for c in range(len(header)):
                _ = 0
                while header[c] != self.header[_]:
                    _ += 1
                self.index.append(_)

            text = self.strand[strand1].bases

            if atindex > end:
                highlighted = str(f"{text[:atindex-1]}[{text[atindex-1:]}]")
            if atindex == end:
                highlighted = str(
                    f"{text[:atindex-1]}[{text[atindex-1:end]}]{text[end:]}"
                )
            else:
                highlighted = str(
                    f"{text[:atindex-1]}[{text[atindex-1:end-1]}]{text[end-1:]}"
                )

            counter = 0
            tools = []
            flag = False

            for counter in range(len(length)):
                if flag:
                    if "]" in highlighted[: length[counter] + 2]:
                        tools.append(highlighted[: length[counter] + 2])
                        highlighted = highlighted[length[counter] + 2 :]
                        flag = False

                    else:
                        tools.append(highlighted[: length[counter] + 1])
                        highlighted = highlighted[length[counter] + 1 :]

                if "[" in highlighted[: length[counter] + 1]:
                    extra = highlighted[: length[counter] + 1]
                    if "]" in highlighted[: length[counter] + 2]:
                        tools.append(highlighted[: length[counter] + 2])
                        highlighted = highlighted[length[counter] + 2 :]
                        flag = False

                    if "[" == extra[-1]:
                        tools.append(highlighted[: length[counter]])
                        highlighted = highlighted[length[counter] :]
                        flag = True

                    else:
                        tools.append(highlighted[: length[counter] + 1])
                        highlighted = highlighted[length[counter] + 1 :]
                        flag = True

                else:
                    tools.append(highlighted[: length[counter]])
                    highlighted = highlighted[length[counter] :]

            if highlighted == "]":
                tools[-1] = tools[-1] + "]"
            while tools[-1] == "":
                tools.pop()

            self.highlighted = self.highlighted + tools
            tools = self.highlighted + tools

            counter = len(self.highlighted) - 1
            while counter >= 0:
                if "[" in self.highlighted[counter] or "]" in self.highlighted[counter]:
                    counter -= 1
                    continue
                tools.pop(counter)
                self.index.pop(counter)
                self.highlighted.pop(counter)
                counter -= 1

            counter = len(self.highlighted) - 1
            while counter >= 0:
                if self.highlighted.count(self.highlighted[counter]) > 1:
                    indices = [
                        i
                        for i, x in enumerate(self.highlighted)
                        if x == self.highlighted[counter]
                    ]
                    while len(indices) > 1:
                        self.highlighted.pop(indices[-1])
                        self.index.pop(indices[-1])
                        indices.pop(-1)
                        counter -= 1
                else:
                    counter -= 1

            for i in range(len(self.thief.columns)):
                if self.thief.loc["energy"][i] == temp["energy"][a]:
                    self.shape.append(temp["shape"][a])
                    text = self.thief.columns.values.tolist()
                    self.thief = self.thief.drop(labels=[text[i]], axis=1)
                    break


def gc_content(gc, self):
    """ "
    This function is used to measure the occurences of Gs or Cs in a given strand and also helps delimit and call specific region identifiers.

    Args:
            (from main class)
            self.strand : Structure class used to define the strand sequences
            self.population_size : Length of a given region provided by the user in user_input()
            self.header : Region identifiers provided by the user in user_input()
            gc : an empty list to store the gc content measured by this function

    Returns:
            data : a dictionary containing region identifiers as keys and their respective sequences as items
            gc : gc content measured throughout each sequence(# of Gs and # of Cs)
    """
    # print(self.parent)
    if self.parent == None:
        gc.clear()
        gci = 0
        temp = []
        for i in range(len(self.strand)):
            temp.append(self.strand[i].bases)

        full_sequence = "".join(temp)

        counter = 0
        temp.clear()
        index = 0
        while counter < len(full_sequence):
            temp.append(full_sequence[counter : counter + self.population_size[index]])
            counter += self.population_size[index]
            index += 1

        for j in range(len(temp)):
            for k in temp[j]:
                if k == "C" or k == "G":
                    gci += 1
            gc.append(gci)
            gci = 0

        data = dict()
        # print(len(self.population_size), len(temp))
        for m in range(len(self.population_size)):
            data[self.header[m]] = temp[m]

    else:
        gc.clear()
        gci = 0
        temp = []
        for i in range(len(self.parent.strand)):
            temp.append(self.parent.strand[i].bases)

        full_sequence = "".join(temp)

        counter = 0
        temp.clear()
        index = 0
        while counter < len(full_sequence):
            temp.append(
                full_sequence[counter : counter + self.parent.population_size[index]]
            )
            counter += self.parent.population_size[index]
            index += 1

        for j in range(len(temp)):
            for k in temp[j]:
                if k == "C" or k == "G":
                    gci += 1
            gc.append(gci)
            gci = 0

        data = dict()
        # print(len(self.parent.population_size), len(temp), len(self.parent.header))
        for m in range(len(self.parent.population_size)):
            data[self.parent.header[m]] = temp[m]

    garbage.collect()
    return data


class Paint_structure(QMainWindow):
    def _init_(self, main):
        super().__init__()
        self.parent = main

    def paint(self, x, y):
        self.setWindowTitle("Structure viewer")
        self._zoom = 0
        self.tabs = QTabWidget()
        self.tabs.setFixedSize(1200, 700)
        self.main_layout = QVBoxLayout()

        cwd = os.getcwd()
        test = os.listdir(cwd)

        for item in test:
            if item.endswith(f"{y}_{x}.png"):
                e = f"{y}_{x}.png"

                label = QLabel(
                    f"This is the interaction between {self.parent().regions[y]} (from 1 to {len(self.parent().strand[y].bases) + 1}) and {self.parent().regions[x]} (from {len(self.parent().strand[y].bases) + 1} to {len(self.parent().strand[y].bases) + len(self.parent().strand[x].bases)})({e})"
                )

                pixmap = QPixmap(e)
                item = QGraphicsPixmapItem(pixmap)
                scene = QGraphicsScene(self)
                scene.addItem(item)
                self.view = QGraphicsView()
                self.view.setScene(scene)

                tab = QWidget()
                tab.layout = QVBoxLayout()
                tab.layout.addWidget(label)
                tab.layout.addWidget(self.view)
                tab.setLayout(tab.layout)
                self.tabs.addTab(
                    tab, f"{self.parent().regions[y]} and {self.parent().regions[x]}"
                )

                self.setCentralWidget(self.tabs)

                self.setFixedSize(1200, 700)
                self.move(100, 100)
                self.show()
                break

    def keyPressEvent(self, event):
        if event.key() == 43:
            self._zoom = 1.25
            self.view.scale(self._zoom, self._zoom)
        elif event.key() == 45:
            self._zoom = 0.8
            self.view.scale(self._zoom, self._zoom)
        elif event.key() == 72:
            self.view.fitInView(self.view.sceneRect(), Qt.KeepAspectRatio)


class DNA_origami(QWidget):
    # This is the main class responsable for the main widget and the different components.

    def __init__(self):
        # View initializer
        super().__init__()
        self.header = []
        self.check = []
        self.fixed_regions = {}
        self.label = []
        self.regions = []
        self.population_size = []
        self.index = []
        self.input_sequence = {}
        self.field = []
        self.raw_structure = None
        self.automate = False
        self.iterator = 0
        self.file_counter = 0
        self.thief = None
        self.strand = []
        self.shape = []
        self.highlighted = []
        self.energy = []
        self.max = 0
        self.parent = None
        self.initUI()

    def initUI(self):
        # This function creates the main components in the main widget, renders them and binds callback functions to components.
        # Set some main window's properties
        self.btn1 = QPushButton("Specify Structure")
        self.btn1.setFixedSize(150, 70)
        self.btn1.setStyleSheet(
            "QPushButton {background-color: #007bff; color: white; border: 5px solid transparent; border-radius: 5px}"
        )
        self.btn1.clicked.connect(self.user_input)
        self.btn2 = QPushButton("Maximum Energy")
        self.btn2.setFixedSize(150, 70)
        self.btn2.setStyleSheet(
            "QPushButton {background-color: #007bff; color: white; border: 5px solid transparent; border-radius: 5px}"
        )
        self.btn2.clicked.connect(self.user_max)
        self.btn3 = QPushButton("Recalculate" + "\n" + "energy")
        self.btn3.setFixedSize(150, 70)
        self.btn3.setIcon(QIcon("assets/sync-solid.svg"))
        self.btn3.setIconSize(QSize(30, 30))
        self.btn3.setStyleSheet(
            "QPushButton {background-color: #007bff; color: white; border: 5px solid transparent; border-radius: 5px}"
        )
        self.btn3.clicked.connect(self.calculate)
        self.btn4 = QPushButton("View/Save" + "\n" + "configuration")
        self.btn4.setFixedSize(150, 70)
        self.btn4.setIcon(QIcon("assets/download-solid.svg"))
        self.btn4.setIconSize(QSize(30, 30))
        self.btn4.setStyleSheet(
            "QPushButton {background-color: #007bff; color: white; border: 5px solid transparent; border-radius: 5px}"
        )
        self.btn4.clicked.connect(self.export)
        self.btn5 = QPushButton("Load" + "\n" + "configuration")
        self.btn5.setFixedSize(150, 70)
        self.btn5.setIcon(QIcon("assets/upload-solid.svg"))
        self.btn5.setIconSize(QSize(30, 30))
        self.btn5.setStyleSheet(
            "QPushButton {background-color: #007bff; color: white; border: 5px solid transparent; border-radius: 5px}"
        )
        self.btn7 = QPushButton()
        self.btn7.setFixedSize(40, 40)
        self.btn7.setToolTip("Export data")
        self.btn7.setIcon(QIcon("assets/save-solid.svg"))
        self.btn7.setIconSize(QSize(20, 20))
        self.btn7.setStyleSheet(
            "QPushButton {border: 5px solid transparent; border-radius: 5px}"
        )
        self.btn8 = QPushButton("Pre-optimize")
        self.btn8.setFixedSize(150, 70)
        self.btn8.setIconSize(QSize(30, 30))
        self.btn8.setStyleSheet(
            "QPushButton {background-color: #007bff; color: white; border: 5px solid transparent; border-radius: 5px}"
        )
        self.btn9 = QPushButton("Calculate" + "\n" + " + " + "\n" + "Pre-optimize")
        self.btn9.setFixedSize(150, 70)
        self.btn9.setIconSize(QSize(30, 30))
        self.btn9.setStyleSheet(
            "QPushButton {background-color: #007bff; color: white; border: 5px solid transparent; border-radius: 5px}"
        )
        self.btn5.clicked.connect(self.load)
        self.btn7.clicked.connect(self.output_data)
        self.btn8.clicked.connect(self.randomize_strand)
        self.btn9.clicked.connect(self.loop_calculation)
        self.tabs = QTabWidget()
        self.tabs.setFixedSize(500, 570)

        self.canvas = FigureCanvas(
            plt.Figure(dpi=100, tight_layout=True, frameon=False)
        )
        self.canvas.setFixedSize(490, 450)
        toolbar = NavigationToolbar2QT(self.canvas, self)
        horizontal_toolbar = QHBoxLayout()
        horizontal_toolbar.addWidget(self.btn7)
        horizontal_toolbar.addWidget(toolbar)
        self.tab1 = QWidget()
        self.tabs.addTab(self.tab1, "Energy")

        self.tab1.layout = QVBoxLayout()
        self.tab1.layout.addWidget(self.canvas)
        self.tab1.layout.addLayout(horizontal_toolbar)
        self.tab1.setLayout(self.tab1.layout)
        self.canvas.mpl_connect("button_press_event", self.render_structure)

        main_layout = QHBoxLayout()

        temp_bottom = QHBoxLayout()
        self.top_layout = QHBoxLayout()
        self.top_layout.setSpacing(150)
        self.center_layout = QFormLayout()
        self.bottom_layout = QHBoxLayout()
        temp_vertical = QVBoxLayout()
        self.bottom_layout.setSpacing(30)
        self.right_side = QVBoxLayout()
        self.left_side = QVBoxLayout()
        self.top_layout.addWidget(self.btn1)
        self.top_layout.addWidget(self.btn2)
        temp_vertical.addWidget(self.btn8)
        temp_vertical.addWidget(self.btn9)
        self.bottom_layout.addLayout(temp_vertical)
        self.bottom_layout.addWidget(self.btn3)
        self.bottom_layout.addWidget(self.btn4)
        self.bottom_layout.addWidget(self.btn5)
        self.left_side.addLayout(self.top_layout)

        self.scroll_1 = QScrollArea()
        self.scroll_1.setStyleSheet("border: none")
        self.mygroup_1 = QGroupBox()
        self.mygroup_1.setStyleSheet("border: none")
        self.mygroup_1.setLayout(self.center_layout)
        self.scroll_1.setWidget(self.mygroup_1)
        self.scroll_1.setWidgetResizable(True)
        self.left_side.addWidget(self.scroll_1)
        self.left_side.addLayout(self.bottom_layout)

        self.right_side.addWidget(self.tabs)
        self.progress = QProgressBar()
        self.progress.setAlignment(Qt.AlignHCenter)
        self.progress.setFixedSize(350, 30)
        temp_bottom.addWidget(self.progress)

        self.right_side.addLayout(temp_bottom)
        main_layout.addLayout(self.left_side)
        main_layout.addLayout(self.right_side)
        self.setLayout(main_layout)
        self.setWindowTitle("DNA origami")
        self.move(100, 100)
        self.setFixedSize(1200, 700)
        self.show()

    def user_input(self):
        """
        This function is used to take the user input and send it to render_form() to create the form layout for the input provided.

        Args:
                self.header : Region identifiers provided by the user in user_input()
                self.label : A text box containing a region identifier with its respective gc content
                self.check : A checkbox to manipulate the status of region identifiers (fixed or modifiable)
                self.field : Fields generated by the user in user_input()
                self.regions : Regions by the user in user_input()
                self.fixed_regions : A dictionary containing the fixed regions as defined by the user in the form of region identifier as key
                                        and associated sequence as item
                self.strand : Structure class used to define the strand sequences
                self.energy : Matrix used to store the energy loss returned by the Mfold program
                self.population_size : Length of a given region provided by the user in user_input()
                self.input_sequence : The randomnly generated sequences from the class Sequence and save for prospective usages
                self.raw_structure : The raw input as inputed by the user from user_input()
        """
        input = QInputDialog()
        text, ok = input.getText(
            self,
            "Input Shape",
            "Enter your desired shape (for example: a25 B25, b25 C25, c25 D25, d25 A25)",
        )
        if ok:
            if str(text) != "":
                while self.center_layout.count() != 0:
                    self.center_layout.removeRow(0)

                self.header.clear()
                self.field.clear()
                self.label.clear()
                self.check.clear()
                self.regions.clear()
                self.fixed_regions = {}
                self.strand.clear()
                self.energy.clear()
                self.population_size.clear()
                self.input_sequence = {}

                try:
                    self.raw_structure = text
                    self.render_form()
                except (KeyError, IndexError):
                    error = QErrorMessage(self)
                    error.showMessage("Please follow the given structure format")

    def render_form(self):
        """
        This function is used to process the input provided by the user and render a form layout containing a header as the region identifier,
        its GC content, a random sequence generated for this region identifier and a checkbox if this region can be modified.

        Args:
                input : User input provided as a string to this function from user_input()
                (from main class)
                self.regions : Regions by the user in user_input()
                self.raw_structure : The raw input as inputed by the user from user_input()
                self.header : Region identifiers provided by the user in user_input()
                self.population_size : Length of a given region provided by the user in user_input()
                self.fixed_regions : A dictionary containing the fixed regions as defined by the user in the form of region identifier as key
                                        and associated sequence as item
                self.strand : Structure class used to define the strand sequences
                self.input_sequences : The randomnly generated sequences from the class Sequence and save for prospective usages
                self.label : A text box containing a region identifier with its respective gc content
                self.check : A checkbox to manipulate the status of region identifiers (fixed or modifiable)
                self.field : Fields generated by the user in user_input()
        """
        self.regions.clear()
        self.header.clear()
        gc = []

        structure = [
            [
                Region(
                    re.findall("\D+", region)[0],
                    int(re.findall("\d+", region)[0]),
                )
                for region in strand
            ]
            for strand in [
                strand.strip().split() for strand in self.raw_structure.split(",")
            ]
        ]

        temp = []

        for i in range(len(structure)):
            region = ""
            for j in range(len(structure[i])):
                self.header.append(str(structure[i][j]).split("'")[1])
                temp.append(int(str(structure[i][j]).split(" ")[1].split(")")[0]))
            for j in range(len(structure[i])):
                if j == 0:
                    region += str(structure[i][j]).split("'")[1]
                else:
                    region += "--" + str(structure[i][j]).split("'")[1]
            self.regions.append(region)
        self.population_size = temp

        sequence = Sequence.random_sequence(structure, self.fixed_regions)
        self.strand = [
            sequence.build_strand(test) for test in sequence.strand_structures
        ]

        if self.input_sequence != {}:
            bases = ""
            for i in range(len(self.regions)):
                for _ in self.regions[i].split("--"):
                    bases += str(self.input_sequence[_])

                self.strand[i].bases = bases
                bases = ""

        data = gc_content(gc, self)

        for i in range(len(self.header)):
            self.label.append(
                QLabel(self.header[i] + " (GC count: " + str(gc[i]) + ") ")
            )
            self.label[i].setFixedSize(170, 30)
            self.check.append(QCheckBox("Set as fixed region"))

            if self.header[i].islower():
                if self.header[i] in self.fixed_regions.keys():
                    self.field.append(QLineEdit(data[self.header[i]]))
                    self.check[i].setChecked(True)
                    self.field[i].setEnabled(False)
                    self.field[i].setMaxLength(self.population_size[i])
                    regex = QRegExp("[ACGT]+")
                    validator = QRegExpValidator(regex)
                    self.field[i].setValidator(validator)
                    horizontal = QHBoxLayout()
                    horizontal.addWidget(self.field[i])
                    horizontal.addWidget(self.check[i])
                else:
                    self.field.append(QLineEdit(data[self.header[i]]))
                    self.field[i].setMaxLength(self.population_size[i])
                    regex = QRegExp("[ACGT]+")
                    validator = QRegExpValidator(regex)
                    self.field[i].setValidator(validator)
                    horizontal = QHBoxLayout()
                    horizontal.addWidget(self.field[i])
                    horizontal.addWidget(self.check[i])
            else:
                self.field.append(QLineEdit(data[self.header[i]]))
                self.field[i].setEnabled(False)
                horizontal = QHBoxLayout()
                horizontal.addWidget(self.field[i])

            self.field[i].setFixedSize(300, 30)
            self.center_layout.addRow(self.label[i], horizontal)
            self.field[i].editingFinished.connect(self.strand_update)
            self.check[i].stateChanged.connect(self.set_fixed)

    def strand_update(self):
        """
        This function is used to fire after the user updates a sequence in the user fields and performs the modifications in the strand variable.

        Args:
                self.field : Fields generated by the user in user_input()
                self.population_size : Length of a given region provided by the user in user_input()
                self.header : Region identifiers provided by the user in user_input()
                self.strand : Structure class used to define the strand sequences
                self.label : A text box containing a region identifier with its respective gc content
        """
        for i in range(len(self.field)):
            self.field[i].setStyleSheet("background-color: white")
            if self.field[i].isModified():
                text = self.field[i].text()

                if len(text) != self.population_size[i]:
                    if not self.automate:
                        error = QErrorMessage(self)
                        error.showMessage(
                            "Please respect the defined structure's length"
                        )
                        data = gc_content([], self)
                        self.field[i].setText(data[self.header[i]])

                else:
                    maj = text.upper()
                    maj = maj[::-1]

                    temp = (
                        maj.replace("A", "temp").replace("T", "A").replace("temp", "T")
                    )
                    temp = (
                        temp.replace("C", "temp").replace("G", "C").replace("temp", "G")
                    )

                    data = gc_content([], self)
                    if self.header[i].islower():
                        for j in range(len(self.strand)):
                            if data[self.header[i].lower()] in self.strand[j].bases:
                                old = self.strand[j].bases
                                update = old.replace(data[self.header[i]], text)
                                self.strand[j].bases = update

                            if self.header[i].upper() in data.keys():
                                if data[self.header[i].upper()] in self.strand[j].bases:
                                    old = self.strand[j].bases
                                    update = old.replace(
                                        data[self.header[i].upper()], temp
                                    )
                                    self.strand[j].bases = update

                    if self.header[i].isupper():
                        for j in range(len(self.strand)):
                            if data[self.header[i].lower()] in self.strand[j].bases:
                                old = self.strand[j].bases
                                update = old.replace(data[self.header[i].lower()], temp)
                                self.strand[j].bases = update

                            if self.header[i].upper() in data.keys():
                                if data[self.header[i].upper()] in self.strand[j].bases:
                                    old = self.strand[j].bases
                                    update = old.replace(
                                        data[self.header[i].upper()], text
                                    )
                                    self.strand[j].bases = update

                    gc = []
                    new = gc_content(gc, self)
                    self.label[i].setText(
                        self.header[i] + " (GC count: " + str(gc[i]) + ") "
                    )

                    for a, b in new.items():
                        if b != data[a]:
                            for c in range(len(self.header)):
                                if a == self.header[c]:
                                    self.label[c].setText(
                                        a + " (GC count: " + str(gc[i]) + ") "
                                    )
                                    self.field[c].setText(b)

    def set_fixed(self):
        """
        This function is used to define fixed regions by the user and store them in a dictionary variable.

        Args:
                self.fixed_regions : A dictionary containing the fixed regions as defined by the user in the form of region identifier as key
                                        and associated sequence as item
                self.check : A checkbox to manipulate the status of region identifiers (fixed or modifiable)
                self.field : Fields generated by the user in user_input()
                self.header : Region identifiers provided by the user in user_input()
        """
        self.fixed_regions = {}
        for i in range(len(self.check)):
            if self.check[i].isChecked() == True:
                self.field[i].setEnabled(False)
                data = gc_content([], self)
                self.fixed_regions[self.header[i]] = data[self.header[i]]
                self.fixed_regions[self.header[i].upper()] = data[
                    self.header[i].upper()
                ]
            if self.check[i].isChecked() == False:
                if self.header[i].islower():
                    self.field[i].setEnabled(True)

    def user_max(self):
        """
        This function is used to specify the max value that is shown in the heatmap and is defined by the user.

        Args:
                self.max : A value set by the user to rescale the heatmap max value
                self.energy : Matrix used to store the energy loss returned by the Mfold program
        """
        energy = QInputDialog()
        text, ok = energy.getText(
            self, "Maximum Energy", "Enter the maximum energy value for this shape:"
        )
        if ok:
            if str(text) != "":
                if float(text) >= 0:
                    self.max = float(text)
                    if self.energy != 0:
                        minimum = []
                        for i in range(len(self.energy)):
                            temp = min(self.energy[i])
                            minimum.append(temp)

                        if self.max > min(minimum):
                            self.update()
                        else:
                            error = QErrorMessage(self)
                            error.showMessage(
                                "This value is less than the minimum possible value"
                            )

                else:
                    error = QErrorMessage(self)
                    error.showMessage("Please give a positive value")

    def export(self):
        """
        This function is used to export the current configuration and respects the upload format required by the widget.

        Args:
                self.regions : Regions provided by the user in user_input()
                self.strand : Structure class used to define the strand sequences
                self.header : Region identifiers provided by the user in user_input()
                self.input_sequence : The randomnly generated sequences from the class Sequence and save for prospective usages
                self.raw_structure : The raw input as inputed by the user from user_input()
                self.fixed_regions : A dictionary containing the fixed regions as defined by the user in the form of region identifier as key
                                        and associated sequence as item
                self.energy : Matrix used to store the energy loss returned by the Mfold program
        """
        message = {}
        data = gc_content([], self)
        ticks = self.regions
        for i in range(len(self.strand)):
            j = len(ticks[i].split("--"))
            for _ in range(j):
                if _ == 0:
                    message[ticks[i]] = str(data[ticks[i].split("--")[_]])
                else:
                    message[ticks[i]] += "--" + str(data[ticks[i].split("--")[_]])

        window = QMessageBox()
        reply = window.question(
            self,
            "Current strand/save current configuration",
            "The current strand configuration is "
            + "\n"
            + str(json.dumps(message, indent=4))
            + "\n"
            + "Do you want to save this configuration ?",
            QMessageBox.Yes | QMessageBox.No,
            0,
        )
        QMessageBox.adjustSize(self)

        if reply == QMessageBox.Yes:
            filename = QFileDialog.getSaveFileName(
                self,
                "Save configuration",
                "/bureau/dna-origami",
            )
            params = {}

            data = gc_content([], self)
            for i in range(len(self.header)):
                self.input_sequence[self.header[i]] = data[self.header[i]]

            params["raw_structure"] = self.raw_structure
            params["mfold_command"] = "./mfold_quik"
            params["boltzmann_factor"] = 1
            params["fixed_regions"] = self.fixed_regions
            params["input_sequence_definitions"] = self.input_sequence
            params["energy_matrix"] = self.energy
            try:
                with open(str(filename[0] + ".dat"), "w") as configfile:
                    yaml.dump(params, configfile)
            except FileNotFoundError:
                pass
        else:
            pass

    def load(self):
        """
        This function is responsable for the upload fonctionality of the widget and takes in a ".dat" file structured as the export format defined.

        Args:
                self.header : Region identifiers provided by the user in user_input()
                self.field : Fields generated by the user in user_input()
                self.label : A text box containing a region identifier with its respective gc content
                self.check : A checkbox to manipulate the status of region identifiers (fixed or modifiable)
                self.regions :  Regions provided by the user in user_input()
                self.fixed_regions : A dictionary containing the fixed regions as defined by the user in the form of region identifier as key
                                        and associated sequence as item
                self.strand : Structure class used to define the strand sequences
                self.population_size : Length of a given region provided by the user in user_input()
                self.input_sequence : The randomnly generated sequences from the class Sequence and save for prospective usages
                self.raw_structure : The raw input as inputed by the user from user_input()
                self.energy : Matrix used to store the energy loss returned by the Mfold program
        """
        filename = QFileDialog.getOpenFileName(
            self,
            "Save configuration",
            "/bureau/dna-origami",
        )

        if filename[0]:
            if not filename[0].endswith(".dat"):
                error = QErrorMessage(self)
                error.showMessage("Please upload a file with .dat extension")
            else:
                try:
                    while self.center_layout.count() != 0:
                        self.center_layout.removeRow(0)

                    self.header.clear()
                    self.field.clear()
                    self.label.clear()
                    self.check.clear()
                    self.regions.clear()
                    self.fixed_regions = {}
                    self.highlighted.clear()
                    self.index.clear()
                    self.strand.clear()
                    self.population_size.clear()
                    self.input_sequence = {}
                    with open(filename[0], "r") as configfile:
                        params = yaml.load(configfile, Loader=yaml.FullLoader)
                        self.raw_structure = params["raw_structure"]
                        self.energy = params["energy_matrix"]
                        self.fixed_regions = params["fixed_regions"]
                        self.input_sequence = params["input_sequence_definitions"]
                        self.render_form()
                        self.update()

                except FileNotFoundError:
                    pass
                except TypeError:
                    error = QErrorMessage(self)
                    error.showMessage("Please upload a file with .dat extension")
                except KeyError:
                    error = QErrorMessage(self)
                    error.showMessage("Please respect the configuration file format")

    def calculate(self):
        """
        This function creates a thread that is responsable for the call to Mfold program and stores the returned energy in a matrix variable.
        It also makes sure that the gc content boundaries are respected.

        Args:
                self.population_size :
                self.strand : Structure class used to define the strand sequences
                self.header : Region identifiers provided by the user in user_input()
                self.fixed_regions : A dictionary containing the fixed regions as defined by the user in the form of region identifier as key
                                        and associated sequence as item
                self.field : Fields generated by the user in user_input()
                self.thread : a thread variable to control the calculation() function
        """
        self.btn3.setEnabled(False)
        flag_range = False
        index = 0

        gc = []
        _ = {}

        try:
            _ = gc_content(gc, self)
        except IndexError:
            self.calculate()

        for j in range(len(gc)):
            if (gc[j] / self.population_size[j]) >= 0.4 and (
                gc[j] / self.population_size[j]
            ) <= 0.6:
                flag_range = True
            else:
                index = j
                flag_range = False
                break

        if self.strand == []:
            error = QErrorMessage(self)
            error.showMessage("Please provide a structure")
            self.btn3.setEnabled(True)

        elif self.header[index] in self.fixed_regions.keys():
            error = QErrorMessage(self)
            error.showMessage(
                f"Please make sure to uncheck the field {self.header[index].lower()} to allow modifications"
            )
            self.btn3.setEnabled(True)

        else:
            if flag_range:
                for i in range(len(self.field)):
                    self.field[i].setToolTip("")

                self.thread = QThread()
                self.worker = Worker(self)
                self.worker.moveToThread(self.thread)

                self.thread.started.connect(self.worker.run_call)

                self.worker.progression.connect(lambda x: self.progress.setValue(x))
                self.worker.progress.connect(lambda x: self.progress.setFormat(x))

                self.worker.render.connect(self.update)
                self.worker.render.connect(self.highlight_thief)

                self.worker.finished.connect(self.worker.deleteLater)
                self.worker.finished.connect(self.terminate_thread)

                self.progress.setTextVisible(True)

                self.thread.start()

            else:
                if gc[index] / self.population_size[index] < 0.4:
                    if self.header[index].islower():
                        temp = self.field[index].text()
                        previous = temp
                        text = temp.replace("A", "G", 1)
                        if previous == text:
                            text = temp.replace("T", "G", 1)
                        self.field[index].setText(text)
                        self.field[index].setModified(True)
                        self.strand_update()
                        self.calculate()

                    else:
                        for j in range(len(self.header)):
                            if self.header[j].islower():
                                if self.header[j] == self.header[index].lower():
                                    index = j
                        temp = self.field[index].text()
                        previous = temp
                        text = temp.replace("A", "G", 1)
                        if previous == text:
                            text = temp.replace("T", "G", 1)
                        self.field[index].setText(text)
                        self.field[index].setModified(True)
                        self.strand_update()
                        self.calculate()

                if gc[index] / self.population_size[index] > 0.6:
                    if self.header[index].islower():
                        temp = self.field[index].text()
                        previous = temp
                        text = temp.replace("G", "A", 1)
                        if previous == text:
                            text = temp.replace("C", "A", 1)
                        self.field[index].setText(text)
                        self.field[index].setModified(True)
                        self.strand_update()
                        self.calculate()

                    else:
                        for j in range(len(self.header)):
                            if self.header[j].islower():
                                if self.header[j] == self.header[index].lower():
                                    index = j
                        temp = self.field[index].text()
                        previous = temp
                        text = temp.replace("G", "A", 1)
                        if previous == text:
                            text = temp.replace("C", "A", 1)
                        self.field[index].setText(text)
                        self.field[index].setModified(True)
                        self.strand_update()
                        self.calculate()

    def terminate_thread(self):
        self.thread.quit()
        self.thread.wait()
        clear_files()

        self.btn3.setEnabled(True)

    def figure_generation(self):
        for e in self.shape:
            _ = e.replace(".png", "").split("_")
            _ = [int(i) for i in _]
            with open(f"{_[0]}_{_[1]}.ct", "r") as configfile:
                counter = 0
                max = 0
                for line in configfile:
                    if line[7:8] == "d":
                        counter += 1
                        continue
                    if counter > (
                        len(self.strand[_[0]].bases) + len(self.strand[_[1]].bases) + 4
                    ):
                        break
                    if int(line[26:30]) > max:
                        max = int(line[26:30])
                    counter += 1

            with open(f"{_[0]}_{_[1]}.ct", "r") as configfile:
                nodes_list = []
                flag = False
                G = nx.MultiDiGraph()

                labels = {}
                text = ""
                counter = 0
                before = 0
                center = 0
                extra = 0
                for line in configfile:
                    if line[7:8] == "d":
                        exit = int(line[3:6])
                        continue

                    if line[7:8] == "L":
                        continue

                    if int(line[26:30]) == max:
                        flag = True

                    if counter < len(self.strand[_[0]].bases):
                        before = int(line[12:15])
                        center = int(line[3:6])
                        nodes_list.append(before)
                        G.add_edge(before, center, color="black", weight=1)

                        if counter == 1:
                            G.add_node(
                                f"{self.regions[_[0]]}:{counter}",
                                weight=f"{self.regions[_[0]]}:{counter}",
                            )
                            labels[f"{self.regions[_[0]]}:{counter}"] = f"{counter}"
                            G.add_edge(
                                f"{self.regions[_[0]]}:{counter}",
                                0,
                                color="black",
                                weight=1,
                            )

                        if (counter + 1) % 10 == 0:
                            G.add_node(
                                f"{self.regions[_[0]]}:{counter+1}",
                                weight=f"{self.regions[_[0]]}:{counter+1}",
                            )
                            labels[f"{self.regions[_[0]]}:{counter+1}"] = f"{counter+1}"
                            G.add_edge(
                                before,
                                f"{self.regions[_[0]]}:{counter+1}",
                                color="black",
                                weight=1,
                            )

                        if int(line[26:30]) != 0:
                            if flag:
                                if int(line[26:30]) < len(self.strand[_[0]].bases):
                                    extra = int(line[26:30])
                                else:
                                    extra = int(line[26:30]) - 3
                                G.add_edge(before, extra, color="red", weight=1)
                            else:
                                if int(line[26:30]) < len(self.strand[_[0]].bases):
                                    extra = int(line[26:30])
                                else:
                                    extra = int(line[26:30]) - 3
                                if not G.has_edge(extra, before):
                                    G.add_edge(before, extra, color="black", weight=1)
                        else:
                            flag = False

                        counter += 1

                    else:
                        before = int(line[12:15]) - 3
                        center = int(line[3:6]) - 3
                        if int(line[26:30]) < len(self.strand[_[0]].bases):
                            extra = int(line[26:30])
                        else:
                            extra = int(int(line[26:30]) - 4)
                        nodes_list.append(before)

                        if (counter + 1) % 10 == 0:
                            G.add_node(
                                f"{self.regions[_[1]]}:{counter+1}",
                                weight=f"{self.regions[_[1]]}:{counter+1}",
                            )
                            labels[f"{self.regions[_[1]]}:{counter+1}"] = f"{counter+1}"
                            G.add_edge(
                                before,
                                f"{self.regions[_[1]]}:{counter+1}",
                                color="black",
                                weight=1,
                            )

                        if exit == int(line[3:6]):
                            if line[7:8] == "d":
                                nodes_list.append(before)
                                G.add_edge(nodes_list[-1], 0, color="black", weight=1)
                                break
                            else:
                                nodes_list.append(before)
                                G.add_edge(nodes_list[-1], 0, color="black", weight=1)
                                break

                        G.add_edge(before, center, color="black", weight=1)
                        if extra != 0:
                            if flag:
                                G.add_edge(before, extra, color="red", weight=1)
                            else:
                                if not G.has_edge(extra, before):
                                    G.add_edge(before, extra, color="black", weight=1)
                        else:
                            flag = False

                        counter += 1

                text = self.strand[_[0]].bases + self.strand[_[1]].bases

                for _ in range(len(text)):
                    G.add_node(nodes_list[_], weight=text[_])
                    labels[_] = G.nodes[_]["weight"]

                options = {
                    "linewidths": 1,
                    "node_color": "#d3d3d3",
                    "width": 1.5,
                    "arrowstyle": "-",
                    "arrowsize": 4,
                    "node_size": 60,
                    "font_size": 5,
                }

                colors = nx.get_edge_attributes(G, "color").values()
                if len(text) < 100:
                    dist = dict(nx.shortest_path_length(G, weight="10"))
                else:
                    dist = dict(nx.shortest_path_length(G, weight=str(10 * len(text))))

                pos = nx.kamada_kawai_layout(G, dist=dist)
                nx.draw(
                    G,
                    pos=pos,
                    with_labels=True,
                    labels=labels,
                    edge_color=colors,
                    **options,
                )
                
                plt.savefig(e, dpi=300)
                plt.close()
                plt.clf()

        cwd = os.getcwd()
        test = os.listdir(cwd)

        for item in test:
            if item.endswith(f".png"):
                if not item in self.shape:
                    os.remove(os.path.join(cwd, item))

        self.shape.clear()

    def highlight_thief(self):
        """
        This function is used to highlighted the fields responsable for the energy loss binding the different strands.

        Args:
                self.highlighted : a list of sequences that are responsable for the energy loss found in highlight_thief()
                self.index : a list of indices corresponding to the sequences in self.highlighted
        """
        self.highlighted.clear()
        self.index.clear()
        self.shape.clear()
        for i in range(2):
            intermediate_process(self)

        if not self.automate:
            t = threading.Thread(target=self.figure_generation)
            t.start()
            t.join()
            del t
            garbage.collect()

        highlighted = []
        highlighted += self.highlighted
        counter = 0
        index = []
        index += self.index
        for d in index:
            if self.header[d].isupper():
                highlighted[counter] = (
                    highlighted[counter][::-1]
                    .replace("C", "temp")
                    .replace("G", "C")
                    .replace("temp", "G")
                )
                highlighted[counter] = (
                    highlighted[counter]
                    .replace("]", "temp")
                    .replace("[", "]")
                    .replace("temp", "[")
                )
                highlighted[counter] = (
                    highlighted[counter]
                    .replace("A", "temp")
                    .replace("T", "A")
                    .replace("temp", "T")
                )

                for _ in range(len(self.header)):
                    if self.header[_] == self.header[d].lower():
                        break
                index[counter] = _
            counter += 1

        flag = False
        counter = 0
        for d in index:
            if counter > len(index):
                break
            if self.field[d].toolTip() != "":
                text = self.field[d].toolTip()
                self.field[d].setToolTip("")
                if flag:
                    if "]" in highlighted[counter]:
                        self.field[d].setToolTip(
                            text + "\n" + highlighted[counter] + "\n"
                        )
                        self.field[d].setStyleSheet("background-color:red")
                        flag = False
                    else:
                        self.field[d].setToolTip(
                            text + "\n" + highlighted[counter] + "\n"
                        )
                        self.field[d].setStyleSheet("background-color:red")

                if "[" in highlighted[counter]:
                    if "]" in highlighted[counter]:
                        self.field[d].setToolTip(
                            text + "\n" + highlighted[counter] + "\n"
                        )

                        self.field[d].setStyleSheet("background-color:red")
                        flag = False
                    else:
                        self.field[d].setToolTip(
                            text + "\n" + highlighted[counter] + "\n"
                        )

                        self.field[d].setStyleSheet("background-color:red")
                        flag = True

            else:
                if flag:
                    if "]" in highlighted[counter]:
                        self.field[d].setToolTip(highlighted[counter] + "\n")

                        self.field[d].setStyleSheet("background-color:red")
                        flag = False
                    else:
                        self.field[d].setToolTip(highlighted[counter] + "\n")
                        self.field[d].setStyleSheet("background-color:red")

                if "[" in highlighted[counter]:
                    if "]" in highlighted[counter]:
                        self.field[d].setToolTip(highlighted[counter] + "\n")

                        self.field[d].setStyleSheet("background-color:red")
                        flag = False
                    else:
                        self.field[d].setToolTip(highlighted[counter] + "\n")
                        self.field[d].setStyleSheet("background-color:red")
                        flag = True
            counter += 1

    def render_structure(self, event):
        self.new_widget = Paint_structure(self)
        try:
            x = int(np.round(event.xdata))
            y = int(np.round(event.ydata))
            self.new_widget.paint(x, y)
        except TypeError:
            pass

    def randomize_strand(self):
        """
        This function is used to automatically correct the energy loss found in calculate().

        Args:
                self.automate : a flag used to bypass and automate calculate() and loop_calculation()
                self.highlighted : a list of sequences that are responsable for the energy loss found in highlight_thief()
                self.index : a list of indices corresponding to the sequences in self.highlighted
                self.header : Region identifiers provided by the user in user_input()
                self.fixed_regions : A dictionary containing the fixed regions as defined by the user in the form of region identifier as key
                                    and associated sequence as item
                self.field : Fields generated by the user in user_input()
        """
        if not self.automate:
            window = QMessageBox()
            reply = window.question(
                self,
                "Randomization",
                "Do you want to randomize the highlighted sequences (GC content may change) ?",
                QMessageBox.Yes | QMessageBox.No,
                0,
            )

        if self.automate:
            reply = QMessageBox.Yes

        if reply == QMessageBox.Yes:
            if self.highlighted == []:
                error = QErrorMessage(self)
                error.showMessage(
                    "Please run the calculation once again to highlight regions"
                )

            else:
                self.btn8.setEnabled(False)
                self.threading = QThread()
                self.work = Randomizer(self)

                self.work.moveToThread(self.threading)

                self.threading.started.connect(self.work.run_randomization)
                self.work.render.connect(self.render_form)
                self.work.progress.connect(lambda x: self.progress.setFormat(x))
                self.work.progression.connect(lambda x: self.progress.setValue(x))
                self.work.error.connect(self.failed_randomization)
                self.work.highlight.connect(self.update)
                self.work.highlight.connect(self.highlight_thief)
                self.work.finished.connect(self.work.deleteLater)
                self.work.finished.connect(self.enabler_8)

                self.threading.start()

    def enabler_8(self):
        self.threading.quit()
        self.threading.wait()
        self.btn8.setEnabled(True)
        del self.threading
        garbage.collect()

        clear_files()

    def failed_randomization(self):
        if self.btn9.isEnabled():
            error = QErrorMessage(self)
            error.showMessage(
                "Randomization failed, please run the randomization once again"
            )

    def loop_calculation(self):
        """
        This function is used to create iteration in the process of lowering energy loss by create nested threads of calculation and randomization.

        Args:
                self.automate : a flag used to bypass and automate calculate() and loop_calculation()
                self.iterator : a user input that determinates how many iterations are needed
                self.file_counter : a counter for the current iteration in loop_calculation()
                self.progress : a progressbar component used in the main layout
                self.threader : a thread variable to control the loop_calculation() function
                self.field : Fields generated by the user in user_input()
        """
        input = QInputDialog()
        text, ok = input.getText(
            self,
            "Iterations required",
            "Please provide the number of iterations",
        )

        if ok:
            if str(text) != "":
                if int(text) < 0:
                    error = QErrorMessage(self)
                    error.showMessage(
                        f"Please make sure to provide an integer value greater than 0"
                    )

                if self.field == []:
                    error = QErrorMessage(self)
                    error.showMessage(f"Please make sure to provide a structure")

                else:
                    if int(text) != 0:
                        self.iterator = int(text)
                        self.btn9.setEnabled(False)
                        self.highlighted.clear()
                        self.index.clear()
                        self.shape.clear()

                        self.automate = True
                        self.threader = QThread()
                        self.looper = Looper(self)

                        self.looper.moveToThread(self.threader)
                        self.threader.started.connect(self.looper.run_routine)
                        self.looper.progress.connect(
                            lambda x: self.progress.setValue(x)
                        )
                        self.looper.highlight.connect(self.highlight_thief)
                        self.looper.render.connect(self.render_form)
                        self.looper.progression.connect(
                            lambda x: self.progress.setFormat(x)
                        )
                        self.looper.finished.connect(self.looper.deleteLater)
                        self.looper.finished.connect(self.update)
                        self.looper.finished.connect(self.highlight_thief)
                        self.looper.finished.connect(self.enabler_9)

                        self.threader.start()

    def enabler_9(self):
        self.threader.quit()
        self.threader.wait()
        self.btn9.setEnabled(True)
        del self.threader

        garbage.collect()
        clear_files()

    def output_data(self):
        """
        This function is used to export the heatmap energy matrix in a ".csv" file.

        Args:
                self.energy : Matrix used to store the energy loss returned by the Mfold program
                self.regions : Region provided by the user in user_input()
        """

        filename = QFileDialog.getSaveFileName(
            self,
            "Save configuration",
            "/bureau/dna-origami",
        )

        if filename[0]:
            if self.energy == []:
                error = QErrorMessage(self)
                error.showMessage(f"There is no heatmap plot at the moment")
            else:
                try:
                    with open(filename[0] + ".csv", "w") as configfile:
                        columns = self.regions
                        df = pd.DataFrame(self.energy, columns=columns)
                        df.to_csv(configfile, index=False)
                except FileNotFoundError:
                    pass

    def update(self):
        """
        This function is used to generate the heatmap from the energy matrix found in calculate().

        Args:
                self.energy : Matrix used to store the energy loss returned by the Mfold program
                self.regions : Region class used to define the regions defined
                self.canvas : A matplotlib figure class to draw the heatmap from the energy matrix stored in calculate()
                self.max : A value set by the user to rescale the heatmap max value
        """
        if self.energy == []:
            error = QErrorMessage(self)
            error.showMessage(f"The energy matrix is empty!")
        else:
            self.canvas.figure.clear()
            ax = self.canvas.figure.add_subplot(111)
            a = self.energy

            ticks = self.regions

            if self.max != 0:
                try:
                    im = ax.imshow(
                        a, cmap="plasma", interpolation="nearest", vmax=self.max
                    )
                except (TypeError, UnboundLocalError):
                    self.update()

            else:
                try:
                    im = ax.imshow(a, cmap="plasma", interpolation="nearest")
                except (TypeError, UnboundLocalError):
                    self.update()

            ax.set_xticks(np.arange(len(ticks)))
            ax.set_yticks(np.arange(len(ticks)))
            ax.set_xticklabels(ticks, fontsize=10, rotation=20, ha="right")
            ax.set_yticklabels(ticks, fontsize=10)
            ax.set_title("Energy Matrix for region interactions")

            ax.figure.colorbar(
                im,
                label="Kcal/mol",
                orientation="horizontal",
            )

            self.canvas.draw()

    def closeEvent(self, event):
        """
        This function is used to close/exit the main widget.
        """
        exit = QMessageBox()
        reply = exit.question(
            self,
            "Confirmation step",
            "Are you sure to quit?",
            QMessageBox.Yes | QMessageBox.No,
            0,
        )
        if reply == QMessageBox.Yes:
            event.accept()
        else:
            event.ignore()


if __name__ == "__main__":
    window = QApplication(sys.argv)
    view = DNA_origami()
    garbage.enable()
    sys.exit(window.exec_())
