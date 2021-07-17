import sys, yaml, json, random, time, os

# Import QApplication and the required widgets from PyQt5.QtWidgets
from PyQt5.QtCore import *
from PyQt5.QtWidgets import *
from PyQt5.QtGui import *
import matplotlib

matplotlib.use("QT5Agg")
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from mfold_library import Region, Mfold
from genetic import Sequence

import re

# abl20 able25, abre25 abr20, ABRE25 abd36 ABLE25, able25 ABD36 abre25, acl20 acle25, acre25 acr20, ACRE25 acd36 ACLE25, acle25 ACD36 acre25
# able25 ABB36 abre25, ABRE25 abb36 ABLE25, acle25 ACB36 acre25, ACRE25 acb36 ACLE25, adle25 ADB36 adre25, ADRE25 adb36 ADLE25

# To do:
# matplotlib resolution image V
# bug fixes
# font size for axis V
# editable highlight only


class Looper(QObject):
    finished = pyqtSignal()
    progression = pyqtSignal()
    progress = pyqtSignal(int)

    def __init__(self, parent):
        QObject.__init__(self)
        self.parent = parent

    def run_routine(self):
        self.parent.calculate()
        while not self.parent.btn3.isEnabled():
            _ = 0

        self.parent.automate = True
        time.sleep(1)
        self.progression.emit()
        time.sleep(1)
        self.parent.automate = True
        self.finished.emit()
        self.progress.emit(int((self.parent.file_counter / self.parent.iterator) * 100))


class Worker(QObject):
    finished = pyqtSignal()
    progress = pyqtSignal(str)
    progression = pyqtSignal(int)

    def __init__(self, parent):
        QObject.__init__(self)
        self.parent = parent

    def run_call(self):
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

        self.progress.emit("DONE")

        find_thief(self, mfold)
        mfold.clean_all()
        self.finished.emit()


def find_thief(self, mfold):
    data = {}
    flag = False
    for i in range(len(self.parent.strand)):
        for j in range(len(self.parent.strand)):
            data[f"{i}{j}"] = {}

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

            else:
                try:
                    with open(f"{i}_{j}.ct", "r") as configfile:
                        counter = 0
                        max = 0

                        for line in configfile:
                            try:
                                if flag and int(line[25:30]) == 0:
                                    if counter <= len(self.parent.strand[i].bases):
                                        end = counter
                                    else:
                                        end = counter - 3
                                    flag = False
                                if flag:
                                    letter += line[7:8]

                                if float(line[25:30]) > max:
                                    letter = ""
                                    max = float(line[25:30])
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
                except FileNotFoundError:
                    data[f"{i}{j}"]["energy"] = 0
                    data[f"{i}{j}"]["letter"] = None
                    data[f"{i}{j}"]["max"] = 0
                    data[f"{i}{j}"]["index_app"] = 0
                    data[f"{i}{j}"]["index_ct"] = 0
                    data[f"{i}{j}"]["end"] = 0
                    pass
    df = pd.DataFrame.from_dict(data)
    self.parent.thief = df


def intermediate_process(self):
    temp = self.thief.transpose()
    for a in range(len(temp["energy"])):
        header = []
        if temp["energy"][a] == temp["energy"].max():
            end = temp["end"][a]
            atindex = temp["index_app"][a]
            if temp["end"][a] - 1 <= len(
                self.strand[int(self.thief.columns[a][0])].bases
            ):
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

            flag = False
            flagger = False
            counter = 0
            for d in self.index:
                tooltip = ""
                if self.field[d].toolTip() != "" and self.index.count(d) > 1:
                    text = self.field[d].toolTip()
                    if flagger:
                        tooltip = f"{text[:end-1]}]{text[end-1:]}"
                    if end > len(text):
                        tooltip = f"{text[: atindex - 1]}[{text[atindex - 1:]}"
                        flagger = True
                    else:
                        tooltip = f"{text[: atindex - 1]}[{text[atindex - 1: end-1]}]{text[end-1:]}"

                    if flag:
                        if "]" in tooltip:
                            self.field[d].setToolTip(tooltip)
                            self.field[d].setStyleSheet("background-color:red")
                            flag = False
                        else:
                            self.field[d].setToolTip(tooltip)
                            self.field[d].setStyleSheet("background-color:red")

                    if "[" in tooltip:
                        if "]" in tooltip:
                            self.field[d].setToolTip(tooltip)
                            self.field[d].setStyleSheet("background-color:red")
                            flag = False
                        else:
                            self.field[d].setToolTip(tooltip)
                            self.field[d].setStyleSheet("background-color:red")
                            flag = True
                    break

                else:
                    if flag:
                        if "]" in tools[counter]:
                            self.field[d].setToolTip(tools[counter])
                            self.field[d].setStyleSheet("background-color:red")
                            flag = False
                        else:
                            self.field[d].setToolTip(tools[counter])
                            self.field[d].setStyleSheet("background-color:red")

                    if "[" in tools[counter]:
                        if "]" in tools[counter]:
                            self.field[d].setToolTip(tools[counter])
                            self.field[d].setStyleSheet("background-color:red")
                            flag = False
                        else:
                            self.field[d].setToolTip(tools[counter])
                            self.field[d].setStyleSheet("background-color:red")
                            flag = True
                    counter += 1

            for i in range(len(self.thief.columns)):
                if self.thief.loc["energy"][i] == temp["energy"][a]:
                    text = self.thief.columns.values.tolist()
                    self.thief = self.thief.drop(labels=[text[i]], axis=1)
                    break


def gc_content(gc, self):
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

    for i in range(len(temp)):
        for j in temp[i]:
            if j == "C" or j == "G":
                gci += 1
        gc.append(gci)
        gci = 0

    data = {}
    for i in range(len(self.population_size)):
        data[self.header[i]] = temp[i]

    return data


def render_form(self, input):
    self.regions.clear()
    gc = []
    if input != None:
        structure = [
            [
                Region(
                    re.findall("\D+", region)[0],
                    int(re.findall("\d+", region)[0]),
                )
                for region in strand
            ]
            for strand in [strand.strip().split() for strand in input.split(",")]
        ]
    else:
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
    self.strand = [sequence.build_strand(test) for test in sequence.strand_structures]

    if input == None:
        bases = ""
        for i in range(len(self.regions)):
            for _ in self.regions[i].split("--"):
                bases += str(self.input_sequence[_])

            self.strand[i].bases = bases
            bases = ""

    data = gc_content(gc, self)

    for i in range(len(self.header)):
        self.label.append(QLabel(self.header[i] + " (GC count: " + str(gc[i]) + ") "))
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


class DNA_origami(QWidget):
    def __init__(self):
        """View initializer."""
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
        self.highlighted = []
        self.energy = []
        self.max = 0
        self.initUI()

    def initUI(self):

        # Set some main window's properties
        self.btn1 = QPushButton("Specify Structure")
        self.btn1.setFixedSize(150, 70)
        self.btn1.setStyleSheet(
            "QPushButton {background-color: #007bff; color: white; border: 5px solid transparent; border-radius: 5px}"
        )
        self.btn1.clicked.connect(self.userinput)
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
        self.btn6 = QPushButton()
        self.btn6.setFixedSize(40, 40)
        self.btn6.setIcon(QIcon("assets/image-solid.svg"))
        self.btn6.setIconSize(QSize(30, 30))
        self.btn6.setStyleSheet(
            "QPushButton {background-color: #fff; color: white; border: 5px solid transparent; border-radius: 5px}"
        )
        self.btn7 = QPushButton()
        self.btn7.setFixedSize(40, 40)
        self.btn7.setIcon(QIcon("assets/save-solid.svg"))
        self.btn7.setIconSize(QSize(30, 30))
        self.btn7.setStyleSheet(
            "QPushButton {background-color: #fff; color: white; border: 5px solid transparent; border-radius: 5px}"
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
        self.btn6.clicked.connect(self.output_image)
        self.btn7.clicked.connect(self.output_data)
        self.btn8.clicked.connect(self.randomize_strand)
        self.btn9.clicked.connect(self.loop_calculation)
        self.canvas = FigureCanvas(plt.Figure(tight_layout=True, frameon=False))
        self.canvas.setFixedSize(500, 500)
        main_layout = QHBoxLayout()

        temp_top = QHBoxLayout()
        temp_bottom = QHBoxLayout()
        temp_top.addWidget(self.btn6)
        temp_top.addWidget(self.btn7)
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
        self.scroll = QScrollArea()
        self.scroll.setStyleSheet("border: none")
        self.mygroup = QGroupBox()
        self.mygroup.setStyleSheet("border: none")
        self.mygroup.setLayout(self.center_layout)
        self.scroll.setWidget(self.mygroup)
        self.scroll.setWidgetResizable(True)
        self.left_side.addWidget(self.scroll)
        self.left_side.addLayout(self.bottom_layout)

        self.right_side.addLayout(temp_top)
        self.right_side.addWidget(self.canvas)
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

    def userinput(self):
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
                    render_form(self, text)
                except KeyError:
                    error = QErrorMessage(self)
                    error.showMessage("Please follow the given structure format")
                except IndexError:
                    error = QErrorMessage(self)
                    error.showMessage("Please follow the given structure format")

    def strand_update(self):
        for i in range(len(self.field)):
            self.field[i].setStyleSheet("background-color: white")
            if self.field[i].isModified():
                text = self.field[i].text()

                if len(text) != self.population_size[i]:
                    error = QErrorMessage(self)
                    error.showMessage("Please respect the defined structure's length")
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
        filename = QFileDialog.getOpenFileName(
            self,
            "Save configuration",
            "/bureau/dna-origami",
        )
        try:
            while self.center_layout.count() != 0:
                self.center_layout.removeRow(0)

            self.header.clear()
            self.field.clear()
            self.label.clear()
            self.check.clear()
            self.regions.clear()
            self.fixed_regions = {}
            self.strand.clear()
            self.population_size.clear()
            self.input_sequence = {}
            with open(filename[0], "r") as configfile:
                params = yaml.load(configfile, Loader=yaml.FullLoader)
                self.raw_structure = params["raw_structure"]
                self.energy = params["energy_matrix"]
                self.fixed_regions = params["fixed_regions"]
                self.input_sequence = params["input_sequence_definitions"]
                render_form(self, None)
                self.update()
        except yaml.parser.ParserError:
            error = QErrorMessage(self)
            error.showMessage("Please upload a file with .dat extension")
        except FileNotFoundError:
            pass
        except TypeError:
            error = QErrorMessage(self)
            error.showMessage("Please upload a file with .dat extension")
        except KeyError:
            error = QErrorMessage(self)
            error.showMessage("Please respect the configuration file format")

    def calculate(self):
        self.btn3.setEnabled(False)
        flag_range = False
        flag_empty = False
        index = 0

        gc = []
        data = gc_content(gc, self)

        for i in range(len(gc)):
            if (gc[i] / self.population_size[i]) >= 0.4 and (
                gc[i] / self.population_size[i]
            ) <= 0.6:
                flag_range = True
            else:
                index = i
                flag_range = False
                break

        if self.strand == []:
            flag_empty = True

        if flag_empty:
            error = QErrorMessage(self)
            error.showMessage("Please provide a structure")
            self.btn3.setEnabled(True)

        if self.header[index] in self.fixed_regions.keys():
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
                self.thread.start()

                self.worker.progression.connect(lambda x: self.progress.setValue(x))
                self.worker.progress.connect(lambda x: self.progress.setFormat(x))
                self.worker.finished.connect(self.update)
                self.worker.finished.connect(self.highlight_thief)
                self.worker.finished.connect(self.thread.quit)
                self.worker.finished.connect(self.worker.deleteLater)
                self.worker.finished.connect(self.thread.wait)
                self.thread.finished.connect(self.thread.deleteLater)

                self.progress.setTextVisible(True)
                self.progress.setFormat("PROGRESSING")

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

    def highlight_thief(self):
        self.highlighted.clear()
        self.index.clear()
        for i in range(2):
            intermediate_process(self)

        self.btn3.setEnabled(True)

    def randomize_strand(self):
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
                flager = True
                flag = False
                while flager:
                    fields = []
                    for i in range(len(self.highlighted)):
                        text = ""

                        for j in range(len(self.highlighted[i])):
                            if self.highlighted[i][j] == "]":
                                flag = False

                            elif flag:
                                text += random.choice(["G", "T", "C", "A"])

                            elif self.highlighted[i][j] == "[":
                                flag = True

                            else:
                                text += self.highlighted[i][j]

                        fields.append(text)

                    counter = len(self.highlighted) - 1
                    while counter >= 0:
                        if (
                            "[" in self.highlighted[counter]
                            or "]" in self.highlighted[counter]
                        ):
                            counter -= 1
                            continue
                        fields.pop(counter)
                        self.index.pop(counter)
                        self.highlighted.pop(counter)
                        counter -= 1

                    gc = []
                    data = gc_content(gc, self)
                    gci = 0
                    for _ in fields[0]:
                        if _ == "C" or _ == "G":
                            gci += 1
                    if gc[self.index[0]] != gci:
                        if "]" in self.highlighted[0]:
                            if "[" in self.highlighted[0]:
                                continue
                            flag = True
                        continue

                    else:
                        if self.header[self.index[0]] in self.fixed_regions.keys():
                            error = QErrorMessage(self)
                            error.showMessage(
                                f"Please make sure to uncheck the field {self.header[self.index[0]].lower()} to allow modifications"
                            )
                            break

                        if self.index.count(self.index[0]) > 1:
                            indices = [
                                i
                                for i, x in enumerate(self.index)
                                if x == self.index[0]
                            ]
                            indices.pop(0)

                            for i in indices:
                                temp = self.highlighted[i]

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
                                temp = [c for c in fields[0]]
                                for k in new_indices_o:
                                    temp.insert(k, "[")

                                for l in new_indices_c:
                                    temp.insert(l, "]")

                                temp = "".join(temp)
                                self.highlighted[i] = temp

                        for i in range(len(self.header)):
                            if self.header[self.index[0]].islower():
                                if self.header[i] == self.header[self.index[0]].upper():
                                    indices = [
                                        j for j, x in enumerate(self.index) if x == i
                                    ]
                                    for k in indices:
                                        temp = self.highlighted[k]
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
                                            fields[0][::-1]
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
                                        self.highlighted[k] = temp

                            if self.header[self.index[0]].isupper():
                                if self.header[i] == self.header[self.index[0]].lower():
                                    indices = [
                                        j for j, x in enumerate(self.index) if x == i
                                    ]
                                    for k in indices:
                                        temp = self.highlighted[k]

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
                                            fields[0][::-1]
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

                                        self.highlighted[k] = temp

                        self.field[self.index[0]].setText(fields[0])
                        self.field[self.index[0]].setModified(True)
                        self.strand_update()
                        self.highlighted.pop(0)
                        self.index.pop(0)
                        try:
                            if "]" in self.highlighted[0]:
                                if "[" in self.highlighted[0]:
                                    continue
                                flag = True
                        except IndexError:
                            pass

                        if self.index == []:
                            break

    def loop_calculation(self):
        if not self.automate:
            input = QInputDialog()
            text, ok = input.getText(
                self,
                "Iterations required",
                "Please provide the number of iterations",
            )

        if self.automate:
            ok = True
            text = self.iterator - self.file_counter
            if text == 0:
                ok = False
                self.automate = False
                self.progress.setValue(100)
            else:
                try:
                    while not self.threader.isFinished():
                        self.threader.wait()
                except RuntimeError:
                    self.automate = False
                    self.loop_calculation()

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
                        if self.file_counter == 0:
                            self.iterator = int(text)
                        self.btn9.setEnabled(False)
                        self.btn8.setEnabled(False)

                        self.automate = True
                        self.threader = QThread()
                        self.looper = Looper(self)
                        self.looper.moveToThread(self.threader)
                        self.threader.started.connect(self.looper.run_routine)

                        self.threader.start()

                        self.looper.progression.connect(self.randomize_strand)
                        self.looper.progression.connect(self.loop_saver)
                        self.looper.progress.connect(
                            lambda x: self.progress.setValue(x)
                        )
                        self.looper.finished.connect(self.threader.quit)
                        self.looper.finished.connect(self.looper.deleteLater)
                        self.looper.finished.connect(self.threader.wait)
                        self.looper.finished.connect(self.loop_calculation)
                        self.threader.finished.connect(self.threader.deleteLater)

        if not ok:
            self.best_run()
            self.iterator = 0
            self.file_counter = 0
            self.automate = False
            self.btn9.setEnabled(True)
            self.btn8.setEnabled(True)

    def output_image(self):
        filename = QFileDialog.getSaveFileName(
            self,
            "Save configuration",
            "/bureau/dna-origami",
        )
        try:
            self.canvas.print_figure(filename[0] + ".png", dpi=300)
        except FileNotFoundError:
            pass

    def output_data(self):
        if self.energy != []:
            filename = QFileDialog.getSaveFileName(
                self,
                "Save configuration",
                "/bureau/dna-origami",
            )
            try:
                with open(filename[0], "w") as configfile:
                    columns = self.regions
                    df = pd.DataFrame(self.energy, columns=columns)
                    df.to_csv(filename[0], index=False)
            except FileNotFoundError:
                pass

    def update(self):
        if self.energy == []:
            error = QErrorMessage(self)
            error.showMessage(f"The energy matrix is empty!")
        else:
            a = self.energy
            i = self.right_side.count()

            ticks = self.regions

            if i != 0:
                self.canvas.figure.clear()

            ax = self.canvas.figure.subplots()

            if self.max != 0:
                im = ax.imshow(a, cmap="hot_r", interpolation="nearest", vmax=self.max)
                ax.set_xticks(np.arange(len(ticks)))
                ax.set_yticks(np.arange(len(ticks)))
                ax.set_xticklabels(
                    ticks, fontsize=10, fontrotation=20, ha="right"
                )
                ax.set_yticklabels(ticks, fontsize=10)
                ax.set_title("Energy Matrix for the different regions")

            else:
                im = ax.imshow(a, cmap="hot_r", interpolation="nearest")
                ax.set_xticks(np.arange(len(ticks)))
                ax.set_yticks(np.arange(len(ticks)))
                ax.set_xticklabels(
                    ticks, fontsize=10, rotation=20, ha="right"
                )
                ax.set_yticklabels(ticks, fontsize=10)
                ax.set_title("Energy Matrix for region interactions")

            cbar = ax.figure.colorbar(im, label="Kcal/mol", orientation="horizontal")
            self.canvas.draw()

    def closeEvent(self, event):
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

    def loop_saver(self):
        if self.file_counter == 0:
            cwd = os.getcwd()
            test = os.listdir(cwd)
            for item in test:
                if item.endswith(".dat"):
                    os.remove(os.path.join(cwd, item))
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
            with open(str(f"loop-{self.file_counter}" + ".dat"), "w") as configfile:
                yaml.dump(params, configfile)
        except FileNotFoundError:
            pass

        self.file_counter += 1

    def best_run(self):
        counter = []
        length = self.file_counter

        for i in range(length):
            higher = 0
            with open(str(f"loop-{i}" + ".dat"), "r") as configfile:
                params = yaml.load(configfile, Loader=yaml.FullLoader)
                matrix = params["energy_matrix"]
                for j in range(len(matrix)):
                    if max(matrix[j]) > higher:
                        higher = max(matrix[j])
                counter.append(higher)

        temp = min(counter)
        index = counter.index(temp)
        error = QErrorMessage(self)
        error.showMessage(
            f"The best iteration with the least energy is in file loop_{index}!"
        )


if __name__ == "__main__":
    window = QApplication(sys.argv)
    view = DNA_origami()
    sys.exit(window.exec_())
