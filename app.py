from io import StringIO
import sys, yaml, json, time

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

from mfold_library import Region, Mfold, EnergyMatrix, Strand
from genetic import Sequence
import re


# To do:
# Sanity checks
# progress bar for energy calculation  medium V
# scroll bar for main panel high V
# scroll bar for save config popup high V
# export energy matrix data medium V
# export energy matrix image medium V
# flag bases contributing to worst energy low
# dashes to separate regions in strand labels high V


def run_call(self):
    mfold = Mfold(output_folder="./", mfold_command="mfold_quik")
    penalty = 1
    energy_matrix = EnergyMatrix(mfold, self.strand, penalty)
    energy_matrix.create()
    self.energy = energy_matrix.matrix
    for i in range(len(self.energy)):
        for j in range(len(self.energy)):
            self.energy[i][j] = abs(self.energy[i][j])

    mfold.clean_all()


def gc_content(gc, self):
    gc.clear()
    gci = 0
    for i in range(len(self.strand)):
        for _ in self.strand[i].bases[: self.population_size[i]]:
            if _ == "G" or _ == "C":
                gci += 1

        gc.append(gci)
        gci = 0
        for _ in self.strand[i].bases[self.population_size[i] :]:
            if _ == "G" or _ == "C":
                gci += 1

        gc.append(gci)
        gci = 0


def render_form(self, input):
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
        self.header.append(str(structure[i][0]).split("'")[1])
        self.header.append(str(structure[i][1]).split("'")[1])
        temp.append(int(str(structure[i][0]).split(" ")[1].split(")")[0]))

    self.population_size = temp

    sequence = Sequence.random_sequence(structure, self.fixed_regions)
    self.strand = [sequence.build_strand(test) for test in sequence.strand_structures]

    if input == None:
        first = []
        last = []
        for i in self.input_sequence.values():
            first.append(i)
            maj = i[::-1]
            temp = maj.replace("A", "temp").replace("T", "A").replace("temp", "T")
            temp = temp.replace("C", "temp").replace("G", "C").replace("temp", "G")
            last.append(temp)

        for i in range(len(self.input_sequence)):
            try:
                self.strand[i].bases = first[i] + last[i + 1]
            except IndexError:
                self.strand[i].bases = first[i] + last[0]

    gc_content(gc, self)

    counter = 0
    for i in range(len(self.header)):
        if i == 0:
            self.label.append(
                QLabel(self.header[i] + " (GC count: " + str(gc[0]) + ") ")
            )
            self.label[i].setFixedSize(150, 30)
            self.check.append(QCheckBox("Set as fixed region"))
        else:
            self.label.append(
                QLabel(self.header[i] + " (GC count: " + str(gc[i]) + ") ")
            )
            self.label[i].setFixedSize(150, 30)
            self.check.append(QCheckBox("Set as fixed region"))

        if i % 2 == 0:
            if self.header[i] in self.fixed_regions.keys():
                self.field.append(
                    QLineEdit(
                        self.strand[counter].bases[: self.population_size[int(i / 2)]]
                    )
                )
                self.check[i].setChecked(True)
                self.field[i].setEnabled(False)
                self.field[i].setMaxLength(self.population_size[int(i / 2)])
                regex = QRegExp("[ACGT]+")
                validator = QRegExpValidator(regex)
                self.field[i].setValidator(validator)
                horizontal = QHBoxLayout()
                horizontal.addWidget(self.field[i])
                horizontal.addWidget(self.check[i])
            else:
                self.field.append(
                    QLineEdit(
                        self.strand[counter].bases[: self.population_size[int(i / 2)]]
                    )
                )
                self.field[i].setMaxLength(self.population_size[int(i / 2)])
                regex = QRegExp("[ACGT]+")
                validator = QRegExpValidator(regex)
                self.field[i].setValidator(validator)
                horizontal = QHBoxLayout()
                horizontal.addWidget(self.field[i])
                horizontal.addWidget(self.check[i])
        else:
            self.field.append(
                QLineEdit(
                    self.strand[counter].bases[self.population_size[int(i / 2)] :]
                )
            )
            self.field[i].setEnabled(False)
            self.field[i].setMaxLength(self.population_size[int(i / 2)])
            regex = QRegExp("[ACGT]+")
            validator = QRegExpValidator(regex)
            self.field[i].setValidator(validator)
            counter += 1
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
        self.population_size = []
        self.input_sequence = {}
        self.field = []
        self.raw_structure = None
        self.strand = []
        self.energy = 0
        self.max = 0
        self.initUI()

    def initUI(self):

        # Set some main window's properties
        self.btn1 = QPushButton("Specify Structure")
        self.btn1.setFixedSize(150, 70)
        self.btn1.setStyleSheet(
            "QPushButton {background-color: #007bff; color: white; border: 5px solid transparent;}"
        )
        self.btn1.clicked.connect(self.userinput)
        self.btn2 = QPushButton("Maximum Energy")
        self.btn2.setFixedSize(150, 70)
        self.btn2.setStyleSheet(
            "QPushButton {background-color: #007bff; color: white; border: 5px solid transparent;}"
        )
        self.btn2.clicked.connect(self.user_max)
        self.btn3 = QPushButton("Recalculate" + "\n" + "energy")
        self.btn3.setFixedSize(150, 70)
        self.btn3.setIcon(QIcon("sync-solid.svg"))
        self.btn3.setIconSize(QSize(30, 30))
        self.btn3.setStyleSheet(
            "QPushButton {background-color: #007bff; color: white; border: 5px solid transparent;}"
        )
        self.btn3.clicked.connect(self.calculate)
        self.btn4 = QPushButton("View/Save" + "\n" + "configuration")
        self.btn4.setFixedSize(150, 70)
        self.btn4.setIcon(QIcon("download-solid.svg"))
        self.btn4.setIconSize(QSize(30, 30))
        self.btn4.setStyleSheet(
            "QPushButton {background-color: #007bff; color: white; border: 5px solid transparent;}"
        )
        self.btn4.clicked.connect(self.export)
        self.btn5 = QPushButton("Load" + "\n" + "configuration")
        self.btn5.setFixedSize(150, 70)
        self.btn5.setIcon(QIcon("upload-solid.svg"))
        self.btn5.setIconSize(QSize(30, 30))
        self.btn5.setStyleSheet(
            "QPushButton {background-color: #007bff; color: white; border: 5px solid transparent;}"
        )
        self.btn6 = QPushButton("Export plot")
        self.btn6.setFixedSize(150, 70)
        self.btn6.setIcon(QIcon("image-solid.svg"))
        self.btn6.setIconSize(QSize(30, 30))
        self.btn6.setStyleSheet(
            "QPushButton {background-color: #007bff; color: white; border: 5px solid transparent;}"
        )
        self.btn7 = QPushButton("Export data")
        self.btn7.setFixedSize(150, 70)
        self.btn7.setIcon(QIcon("save-solid.svg"))
        self.btn7.setIconSize(QSize(30, 30))
        self.btn7.setStyleSheet(
            "QPushButton {background-color: #007bff; color: white; border: 5px solid transparent;}"
        )
        self.btn5.clicked.connect(self.load)
        self.btn6.clicked.connect(self.output_image)
        self.btn7.clicked.connect(self.output_data)
        self.canvas = FigureCanvas(plt.Figure(figsize=(8, 8)))
        self.canvas.setFixedSize(500, 500)
        main_layout = QHBoxLayout()

        temp_top = QHBoxLayout()
        temp_top.addWidget(self.btn6)
        temp_top.addWidget(self.btn7)
        self.top_layout = QHBoxLayout()
        self.top_layout.setSpacing(150)
        self.center_layout = QFormLayout()
        self.bottom_layout = QHBoxLayout()
        self.bottom_layout.setSpacing(50)
        self.right_side = QVBoxLayout()
        self.left_side = QVBoxLayout()
        self.top_layout.addWidget(self.btn1)
        self.top_layout.addWidget(self.btn2)
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
        # progress.setGeometry(100, 30, 350, 30)
        self.progress.setFixedSize(350, 30)
        self.right_side.addWidget(self.progress)
        # self.right_side.addWidget(self.btn6)
        # self.right_side.addWidget(self.btn7)
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
                self.fixed_regions = {}
                self.strand.clear()
                self.population_size.clear()
                self.header.clear()

                try:
                    render_form(self, text)
                except KeyError:
                    error = QErrorMessage(self)
                    error.showMessage("Please follow the given structure format")
                except IndexError:
                    error = QErrorMessage(self)
                    error.showMessage("Please follow the given structure format")

    def strand_update(self):
        flag = False
        for i in range(len(self.field)):
            if self.field[i].isModified():
                text = self.field[i].text()

                if len(text) != self.population_size[int(i / 2)]:
                    flag = True

                if flag == False:
                    self.field[i].setText(text.upper())

                    last = self.strand[int(i / 2)].bases[
                        self.population_size[int(i / 2)] :
                    ]

                    self.strand[int(i / 2)].bases = text.upper() + last

                    gc = []
                    gc_content(gc, self)
                    self.label[i].setText(
                        self.header[i] + " (GC count: " + str(gc[i]) + ") "
                    )
                    maj = text.upper()
                    maj = maj[::-1]

                    temp = (
                        maj.replace("A", "temp").replace("T", "A").replace("temp", "T")
                    )
                    temp = (
                        temp.replace("C", "temp").replace("G", "C").replace("temp", "G")
                    )

                    if i == 0:
                        first = self.strand[-1].bases[
                            : self.population_size[int(i / 2)]
                        ]
                        self.strand[-1].bases = first + temp
                        gc_content(gc, self)
                        self.field[-1].setText(temp)
                        self.label[-1].setText(
                            self.header[-1] + " (GC count: " + str(gc[-1]) + ") "
                        )

                    else:
                        first = self.strand[int((i / 2) - 1)].bases[
                            : self.population_size[int(i / 2)]
                        ]
                        self.strand[int(i / 2) - 1].bases = first + temp
                        gc_content(gc, self)
                        self.field[i - 1].setText(temp)
                        self.label[i - 1].setText(
                            self.header[i - 1] + " (GC count: " + str(gc[i]) + ") "
                        )

            while flag:
                error = QErrorMessage(self)
                error.showMessage("Please respect the defined structure's length")
                self.field[i].setText(
                    self.strand[int(i / 2)].bases[: self.population_size[int(i / 2)]]
                )
                flag = False

    def set_fixed(self):
        self.fixed_regions = {}
        for i in range(len(self.check)):
            if self.check[i].isChecked() == True:
                self.field[i].setEnabled(False)
                self.fixed_regions[self.header[i]] = self.strand[int(i / 2)].bases[
                    : self.population_size[int(i / 2)]
                ]
            if self.check[i].isChecked() == False:
                if i % 2 == 0:
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
        ticks = [
            self.header[i] + "--" + self.header[i + 1]
            for i in range(0, len(self.header), 2)
        ]
        for i in range(len(self.strand)):
            message[ticks[i]] = str(
                self.strand[i].bases[: self.population_size[i]]
                + "--"
                + self.strand[i].bases[self.population_size[i] :]
            )

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
            for i in range(len(self.header)):
                if i % 2 == 0:
                    self.input_sequence[self.header[i]] = self.strand[int(i / 2)].bases[
                        : self.population_size[int(i / 2)]
                    ]

            params["raw_structure"] = self.raw_structure
            params["mfold_command"] = "./mfold_quik"
            params["boltzmann_factor"] = 1
            params["fixed_regions"] = self.fixed_regions
            params["input_sequence_definitions"] = self.input_sequence
            params["energy_matrix"] = self.energy
            try:
                with open(filename[0], "w") as configfile:
                    yaml.dump(params, configfile)
            except FileNotFoundError:
                pass
        else:
            pass

    def load(self):
        while self.center_layout.count() != 0:
            self.center_layout.removeRow(0)
        self.header.clear()
        self.field.clear()
        self.fixed_regions = {}
        self.check.clear()
        self.label.clear()
        self.strand.clear()

        filename = QFileDialog.getOpenFileName(
            self,
            "Save configuration",
            "/bureau/dna-origami",
        )
        try:
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
        except KeyError:
            error = QErrorMessage(self)
            error.showMessage("Please respect the configuration file format")

    def calculate(self):
        flag_range = False
        self.progress.setFormat("PROCESSING")

        gc = []
        gc_content(gc, self)
        for i in range(len(gc)):
            if (
                gc[i] / self.population_size[int(i / 2)] < 0.4
                or gc[i] / self.population_size[int(i / 2)] > 0.6
            ):
                flag_range = True

        if flag_range:
            error = QErrorMessage(self)
            error.showMessage("Please stay in the 40%-60% range for the gc count")
            flag_range = False

        else:
            self.progress.setTextVisible(True)
            self.progress.setValue(100)
            self.progress.setFormat("PROCESSING")
            self.progress.setAlignment(Qt.AlignCenter)

            run_call(self)
            self.progress.setFormat("DONE")
            self.progress.setAlignment(Qt.AlignCenter)
            try:
                self.update()
            except TypeError:
                error = QErrorMessage(self)
                error.showMessage("Please give a structure")

            self.progress.setValue(0)

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
                    columns = [
                        self.header[i] + self.header[i + 1]
                        for i in range(0, len(self.header), 2)
                    ]
                    df = pd.DataFrame(self.energy, columns=columns)
                    df.to_csv(filename[0], index=False)
            except FileNotFoundError:
                pass

    def update(self):
        a = self.energy
        i = self.right_side.count()

        ticks = [
            self.header[i] + self.header[i + 1] for i in range(0, len(self.header), 2)
        ]

        if i != 0:
            self.canvas.figure.clear()

        ax = self.canvas.figure.subplots()

        if self.max != 0:
            im = ax.imshow(a, cmap="hot", interpolation="nearest", vmax=self.max)
            ax.set_xticks(np.arange(len(ticks)))
            ax.set_yticks(np.arange(len(ticks)))
            ax.set_xticklabels(ticks)
            ax.set_yticklabels(ticks)
            ax.tick_params(top=True, bottom=False, labeltop=True, labelbottom=False)
            ax.set_title("Energy Matrix for the different regions")

        else:
            im = ax.imshow(a, cmap="hot", interpolation="nearest")
            ax.set_xticks(np.arange(len(ticks)))
            ax.set_yticks(np.arange(len(ticks)))
            ax.set_xticklabels(ticks)
            ax.set_yticklabels(ticks)
            ax.tick_params(top=True, bottom=False, labeltop=True, labelbottom=False)
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
