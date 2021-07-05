import sys, yaml, json, random

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


# To do:
# Sanity checks
# GC content flag automated V
# Pre optimize highlighted fields
# sanity checks highlighted fields (a25... )


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

        self.progression.emit(100)

        self.progress.emit("DONE")

        find_thief(self, mfold)
        mfold.clean_all()
        self.finished.emit()
        self.parent.btn3.setEnabled(True)


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
        for _ in temp[i]:
            if _ == "G" or _ == "C":
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
        self.btn5.clicked.connect(self.load)
        self.btn6.clicked.connect(self.output_image)
        self.btn7.clicked.connect(self.output_data)
        self.btn8.clicked.connect(self.randomize_strand)
        self.canvas = FigureCanvas(plt.Figure(figsize=(7, 7)))
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
        self.bottom_layout.setSpacing(30)
        self.right_side = QVBoxLayout()
        self.left_side = QVBoxLayout()
        self.top_layout.addWidget(self.btn1)
        self.top_layout.addWidget(self.btn2)
        self.bottom_layout.addWidget(self.btn8)
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
                self.population_size.clear()
                self.input_sequence = {}

                try:
                    render_form(self, text)
                    self.raw_structure = text
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

        except FileNotFoundError:
            pass
        except TypeError:
            error = QErrorMessage(self)
            error.showMessage("Please respect the configuration file type")
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
                self.worker.finished.connect(self.thread.quit)
                self.worker.finished.connect(self.update)
                self.worker.finished.connect(self.highlight_thief)
                self.worker.finished.connect(self.worker.deleteLater)
                self.thread.finished.connect(self.thread.deleteLater)

                self.progress.setValue(1)
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
        header = []
        temp = self.thief.transpose()
        self.highlighted.clear()
        self.index.clear()

        for i in range(len(temp["energy"])):
            header.clear()
            if temp["energy"][i] == temp["energy"].max():
                end = temp["end"][i]
                atindex = temp["index_app"][i]
                if temp["end"][i] - 1 <= len(
                    self.strand[int(self.thief.columns[i][0])].bases
                ):
                    strand1 = int(self.thief.columns[i][0])
                else:
                    strand1 = int(self.thief.columns[i][1])
                    end = end - len(self.strand[int(self.thief.columns[i][0])].bases)
                    atindex = atindex - len(
                        self.strand[int(self.thief.columns[i][0])].bases
                    )

                for j in range(len(self.regions[strand1].split("--"))):
                    header.append(self.regions[strand1].split("--")[j])

                data = gc_content([], self)
                length = [len(data[header[_]]) for _ in range(len(header))]

                for i in range(len(header)):
                    _ = 0
                    while header[i] != self.header[_]:
                        _ += 1
                    self.index.append(_)

                text = self.strand[strand1].bases
                highlighted = str(
                    f"{text[:atindex-1]}[{text[atindex-1:end-1]}]{text[end-1:]}"
                )

                tools = []
                counter = 0
                flag = False
                end = False
                start = True
                for counter in range(len(length)):
                    if end:
                        tools.append(highlighted[: length[counter] + 2])
                        highlighted = highlighted[length[counter] + 2 :]
                        flag = False

                    if "]" in highlighted[: length[counter] + 2]:
                        tools.append(highlighted[: length[counter] + 2])
                        highlighted = highlighted[length[counter] + 2 :]
                        flag = False
                        end = True

                    if flag:
                        tools.append(highlighted[: length[counter] + 1])
                        highlighted = highlighted[length[counter] + 1 :]

                    if "[" in highlighted[: length[counter] + 1]:
                        if "]" in highlighted[: length[counter] + 2]:
                            tools.append(highlighted[: length[counter] + 2])
                            highlighted = highlighted[length[counter] + 2 :]
                            flag = False
                            start = False
                            end = True
                        else:
                            tools.append(highlighted[: length[counter] + 1])
                            highlighted = highlighted[length[counter] + 1 :]
                            flag = True
                            start = False

                    if start:
                        tools.append(highlighted[: length[counter]])
                        highlighted = highlighted[length[counter] :]

                if highlighted == "]":
                    tools[-1] = tools[-1] + "]"
                while tools[-1] == "":
                    tools.pop()

                self.highlighted.append(tools)

                counter = 0
                for i in self.index:
                    if flag:
                        if "]" in tools[counter]:
                            self.field[i].setToolTip(tools[counter])
                            self.field[i].setStyleSheet("background-color:red")
                            break
                        else:
                            self.field[i].setToolTip(tools[counter])
                            self.field[i].setStyleSheet("background-color:red")
                            counter += 1

                    if "[" in tools[counter]:
                        if "]" in tools[counter]:
                            self.field[i].setToolTip(tools[counter])
                            self.field[i].setStyleSheet("background-color:red")
                            break
                        else:
                            self.field[i].setToolTip(tools[counter])
                            self.field[i].setStyleSheet("background-color:red")
                            flag = True
                            counter += 1

                    else:
                        counter += 1
                break

    def randomize_strand(self):
        window = QMessageBox()
        reply = window.question(
            self,
            "Randomization",
            "Do you want to randomize the highlighted sequences (GC content may change) ?",
            QMessageBox.Yes | QMessageBox.No,
            0,
        )
        if reply == QMessageBox.Yes:
            if self.highlighted == []:
                error = QErrorMessage(self)
                error.showMessage(
                    "Please run the calculation once again to highlight regions"
                )

            else:
                fields = []
                for i in range(len(self.highlighted)):
                    flag = False
                    for j in range(len(self.highlighted[i])):
                        text = ""

                        for k in range(len(self.highlighted[i][j])):
                            if self.highlighted[i][j][k] == "]":
                                flag = False

                            elif flag:
                                text += random.choice("ACGT")

                            elif self.highlighted[i][j][k] == "[":
                                flag = True

                            else:
                                text += self.highlighted[i][j][k]

                        if text != "":
                            fields.append(text)

                # abl20 able25, abre25 abr20, ABRE25 abd36 ABLE25, able25 ABD36 abre25, acl20 acle25, acre25 acr20, ACRE25 acd36 ACLE25, acle25 ACD36 acre25

                counter = 0
                for i in self.index:
                    if self.header[i] in self.fixed_regions.keys():
                        error = QErrorMessage(self)
                        error.showMessage(
                            f"Please make sure to uncheck the field {self.header[i].lower()} to allow modifications"
                        )
                        break
                    else:
                        self.field[i].setText(fields[counter])
                        self.field[i].setModified(True)
                        counter += 1
                        self.strand_update()

                if counter == len(self.index):
                    self.highlighted.clear()

    def output_image(self):
        filename = QFileDialog.getSaveFileName(
            self,
            "Save configuration",
            "/bureau/dna-origami",
        )
        try:
            self.canvas.print_png(filename[0])
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
                im = ax.imshow(a, cmap="hot", interpolation="nearest", vmax=self.max)
                ax.set_xticks(np.arange(len(ticks)))
                ax.set_yticks(np.arange(len(ticks)))
                ax.set_xticklabels(ticks, rotation=20, ha="right")
                ax.set_yticklabels(ticks)
                ax.set_title("Energy Matrix for the different regions")

            else:
                im = ax.imshow(a, cmap="hot", interpolation="nearest")
                ax.set_xticks(np.arange(len(ticks)))
                ax.set_yticks(np.arange(len(ticks)))
                ax.set_xticklabels(ticks, rotation=20, ha="right", fontsize=8)
                ax.set_yticklabels(ticks, fontsize=8)
                ax.set_title("Energy Matrix for the different regions")

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


if __name__ == "__main__":
    window = QApplication(sys.argv)
    view = DNA_origami()
    sys.exit(window.exec_())
