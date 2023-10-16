"""
WIP
Planar Tab of Main Window
- 2D Array Factors
- Elevation Plane and Azimuth Plane Graphing
- 3D Rectangular Coordinate Plotting
"""
from PySide6.QtWidgets import QLabel, QCheckBox, QWidget, QFormLayout, QHBoxLayout, QVBoxLayout, QButtonGroup, QTabWidget, QLineEdit, QSlider, QPushButton, QMessageBox
from PySide6.QtCore import Qt

import matplotlib
import matplotlib.pyplot as plt
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg

import numpy as np


# Widget in the Planar Array tab
class PlanarArrayTab(QWidget):
    def __init__(self):    # Widget Constructor
        super().__init__()
        # Initial Array Parameters
        self.N_elements1 = 2
        self.N_elements2 = 2
        self.axis = "XY"
        self.distance1 = 2
        self.distance2 = 2
        self.phaseshift1 = 0
        self.phaseshift2 = 0
        self.normalize = True

        # Number of Array Elements Axis 1
        self.NumberofElements1 = QLineEdit()
        self.NumberofElements1.setMaximumHeight(30)
        self.NumberofElements1.setMaximumWidth(40)
        self.NumberofElements1.setText("2")
        self.NumberofElements1.returnPressed.connect(self.setElementNum1)
        self.NumberofElements1.returnPressed.connect(self.updateGraphs)

        # Number of Array Elements Axis 2
        self.NumberofElements2 = QLineEdit()
        self.NumberofElements2.setMaximumHeight(30)
        self.NumberofElements2.setMaximumWidth(40)
        self.NumberofElements2.setText("2")
        self.NumberofElements2.returnPressed.connect(self.setElementNum2)
        self.NumberofElements2.returnPressed.connect(self.updateGraphs)

        # Distance Horizontal Layout
        enterNumElementsLayout = QFormLayout()
        enterNumElementsLayout.addRow("Enter Number of Elements (Axis1):", self.NumberofElements1)
        enterNumElementsLayout.addRow("Enter Number of Elements (Axis2):", self.NumberofElements2)
        # =========================================================================================
        self.AxisLabel = QLabel("Choose Axis:")
        # Axis Checkboxes
        self.chooseXYAxis = QCheckBox("XY")
        self.chooseXYAxis.setChecked(True)
        self.chooseXZAxis = QCheckBox("XZ")
        self.chooseYZAxis = QCheckBox("YZ")

        # Grouping Axis Checkboxes
        self.AxisGroup = QButtonGroup()
        self.AxisGroup.addButton(self.chooseXYAxis)
        self.AxisGroup.addButton(self.chooseXZAxis)
        self.AxisGroup.addButton(self.chooseYZAxis)
        self.AxisGroup.buttonClicked.connect(self.setAxis)
        self.AxisGroup.buttonClicked.connect(self.updateGraphs)

        # Axis Horizontal Layout
        chooseAxisLayout = QHBoxLayout()
        chooseAxisLayout.addWidget(self.chooseXYAxis)
        chooseAxisLayout.addWidget(self.chooseXZAxis)
        chooseAxisLayout.addWidget(self.chooseYZAxis)
        # =========================================================================================
        self.AmplitudeLabel = QLabel("Choose Amplitude Excitation:")
        # Amplitude Checkboxes
        self.chooseUniformAmplitude = QCheckBox("Uniform")
        self.chooseUniformAmplitude.setChecked(True)
        self.chooseVariedAmplitude = QCheckBox("Variable Gain")
        self.chooseVariedAmplitude.toggled.connect(self.ampExcitationToggle)  # Toggles Excitation List Layout to input amplitudes

        # Grouping Amplitude Checkboxes
        self.AmplitudeGroup = QButtonGroup()
        self.AmplitudeGroup.addButton(self.chooseUniformAmplitude)
        self.AmplitudeGroup.addButton(self.chooseVariedAmplitude)

        # Variable Amplitude Text Edit
        self.AmplitudeExcitationLabel = QLabel("Excitations:")
        self.UniformAmpLabel = QLabel("Unity")
        self.AmplitudeExcitationLabel.hide()
        # Amplitudes List 2
        self.AmplitudeExcitationList1 = QLineEdit()
        self.AmplitudeExcitationList1.setPlaceholderText("a1, a2, a3, a4, ...")
        self.AmplitudeExcitationList1.hide()
        self.AmplitudeExcitationList1.setMaximumHeight(30)
        self.AmplitudeExcitationList1.setMaximumWidth(100)
        # Amplitudes List 2
        self.AmplitudeExcitationList2 = QLineEdit()
        self.AmplitudeExcitationList2.setPlaceholderText("a1, a2, a3, a4, ...")
        self.AmplitudeExcitationList2.hide()
        self.AmplitudeExcitationList2.setMaximumHeight(30)
        self.AmplitudeExcitationList2.setMaximumWidth(100)

        # Amplitude Text Edit Layout
        self.ampExcitationList = QFormLayout()
        self.ampExcitationList.addRow("Excitations (Axis1):", self.AmplitudeExcitationList1)
        self.ampExcitationList.addRow("Excitations (Axis2):", self.AmplitudeExcitationList2)
        self.ampExcitationList.addRow("All Excitations =", self.UniformAmpLabel)
        self.ampExcitationList.setRowVisible(self.AmplitudeExcitationList1, False)
        self.ampExcitationList.setRowVisible(self.AmplitudeExcitationList2, False)
        self.ampExcitationList.setRowVisible(self.UniformAmpLabel, True)

        # Amplitude Horizontal Layout
        chooseAmplitudeLayout = QHBoxLayout()
        chooseAmplitudeLayout.addWidget(self.chooseUniformAmplitude)
        chooseAmplitudeLayout.addWidget(self.chooseVariedAmplitude)
        # =========================================================================================
        # Distance Text Edit
        self.Distance1Entry = QLineEdit()
        self.Distance1Entry.setMaximumHeight(30)
        self.Distance1Entry.setMaximumWidth(30)
        self.Distance1Entry.setText("2")
        self.Distance1Entry.returnPressed.connect(self.setElementDist)
        self.Distance1Entry.returnPressed.connect(self.updateGraphs)

        self.Distance2Entry = QLineEdit()
        self.Distance2Entry.setMaximumHeight(30)
        self.Distance2Entry.setMaximumWidth(30)
        self.Distance2Entry.setText("2")
        self.Distance2Entry.returnPressed.connect(self.setElementDist)
        self.Distance2Entry.returnPressed.connect(self.updateGraphs)


        # Distance Horizontal Layout
        enterDistanceLayout = QFormLayout()
        enterDistanceLayout.addRow("Enter Element Distance (Axis1): λ/", self.Distance1Entry)
        enterDistanceLayout.addRow("Enter Element Distance (Axis2): λ/", self.Distance2Entry)
        # =========================================================================================
        self.PhaseShiftLabel = QLabel("Select Progressive Phase Shift:")
        self.PhaseValLabel1 = QLabel("β_1 = 0")
        self.PhaseValLabel2 = QLabel("β_2 = 0")
        # Phase Shift Slider 1
        self.PhaseSlider1 = QSlider()
        self.PhaseSlider1.setMinimum(-100)
        self.PhaseSlider1.setMaximum(100)
        self.PhaseSlider1.setTickInterval(1)
        self.PhaseSlider1.setOrientation(Qt.Orientation.Horizontal)
        self.PhaseSlider1.valueChanged.connect(self.updatePhaseValLabel1)
        self.PhaseSlider1.valueChanged.connect(self.updateGraphs)
        # Phase Shift Slider 2
        self.PhaseSlider2 = QSlider()
        self.PhaseSlider2.setMinimum(-100)
        self.PhaseSlider2.setMaximum(100)
        self.PhaseSlider2.setTickInterval(1)
        self.PhaseSlider2.setOrientation(Qt.Orientation.Horizontal)
        self.PhaseSlider2.valueChanged.connect(self.updatePhaseValLabel2)
        self.PhaseSlider2.valueChanged.connect(self.updateGraphs)
        # Phase Shift Sweep Select
        self.PhaseShiftSweep1 = QCheckBox("Sweep Phase 1")
        self.PhaseShiftSweep2 = QCheckBox("Sweep Phase 2")
        # Resets Slider to β = 0
        self.ResetSlider1 = QPushButton("Reset Slider 1")
        self.ResetSlider1.clicked.connect(self.resetSlider1)
        self.ResetSlider2 = QPushButton("Reset Slider 2")
        self.ResetSlider2.clicked.connect(self.resetSlider2)

        # Button and Checkbox Layout
        phasebuttons_layout1 = QHBoxLayout()
        phasebuttons_layout1.addWidget(self.PhaseShiftSweep1)
        phasebuttons_layout1.addWidget(self.PhaseShiftSweep2)
        phasebuttons_layout2 = QHBoxLayout()
        phasebuttons_layout2.addWidget(self.ResetSlider1)
        phasebuttons_layout2.addWidget(self.ResetSlider2)
        phasebuttons_layout = QVBoxLayout()
        phasebuttons_layout.addLayout(phasebuttons_layout1)
        phasebuttons_layout.addLayout(phasebuttons_layout2)
        # =========================================================================================
        # Normalize AF checkbox
        self.NormalizeCheck = QCheckBox("Normalize")
        self.NormalizeCheck.toggled.connect(self.setNormalized)
        self.NormalizeCheck.setChecked(True)
        # =========================================================================================
        # 2D and 3D Graphing
        #self.TwoDimensional_Plt = LinearElevationAzimuthPlot()
        #self.TwoDimensional_Plt.graphAF(self.axis, self.N_elements, self.distance, self.phaseshift, self.normalize)

        #self.ThreeDimensional_Plt = LinearFull3DPlot()
        #self.ThreeDimensional_Plt.graph3DAF(self.axis, self.N_elements, self.distance, self.phaseshift, self.normalize)

        # Tabs for Different Graphs
        self.graph_tab_menu = QTabWidget()
        #self.graph_tab_menu.addTab(, "2D Plots")
        #self.graph_tab_menu.addTab(, "3D Plot")
        # =========================================================================================
        # Configuration Layout
        config_layout = QVBoxLayout()
        config_layout.addWidget(self.AxisLabel)
        config_layout.addLayout(chooseAxisLayout)
        config_layout.addLayout(enterNumElementsLayout)
        config_layout.addWidget(self.AmplitudeLabel)
        config_layout.addLayout(chooseAmplitudeLayout)
        config_layout.addLayout(self.ampExcitationList)
        config_layout.addLayout(enterDistanceLayout)
        config_layout.addWidget(self.PhaseShiftLabel)
        config_layout.addWidget(self.PhaseValLabel1)
        config_layout.addWidget(self.PhaseSlider1)
        config_layout.addWidget(self.PhaseValLabel2)
        config_layout.addWidget(self.PhaseSlider2)
        config_layout.addLayout(phasebuttons_layout)
        config_layout.addWidget(self.NormalizeCheck)

        # Main Layout Combining Config and Graphs
        main_layout = QHBoxLayout()
        main_layout.addLayout(config_layout)
        main_layout.addWidget(self.graph_tab_menu)
        main_layout.setStretch(1, 1.1)

        self.setLayout(main_layout)

    def setElementNum1(self):
        val = self.NumberofElements1.text()
        try:
            val = int(val)
        except ValueError:
            pass
        if isinstance(val, int) and val > 0:
            self.N_elements1 = val
            return
        else:
            N_elements_alert = QMessageBox()
            N_elements_alert.setWindowTitle("Bad Input")
            N_elements_alert.setText("Entry is not an integer, or <= 0")
            N_elements_alert.setInformativeText("Retry")
            N_elements_alert.setIcon(QMessageBox.Icon.Critical)
            N_elements_alert.exec()

    def setElementNum2(self):
        val = self.NumberofElements2.text()
        try:
            val = int(val)
        except ValueError:
            pass
        if isinstance(val, int) and val > 0:
            self.N_elements2 = val
            return
        else:
            N_elements_alert = QMessageBox()
            N_elements_alert.setWindowTitle("Bad Input")
            N_elements_alert.setText("Entry is not an integer, or <= 0")
            N_elements_alert.setInformativeText("Retry")
            N_elements_alert.setIcon(QMessageBox.Icon.Critical)
            N_elements_alert.exec()

    def setElementDist(self):
        val1 = self.Distance1Entry.text()
        val2 = self.Distance2Entry.text()
        try:
            val1 = float(val1)
            val2 = float(val2)
        except ValueError:
            pass
        if (isinstance(val1, float) and val1 > 0) and (isinstance(val2, float) and val2 > 0):
            self.distance1 = val1
            self.distance2 = val2
            return
        else:
            N_elements_alert = QMessageBox()
            N_elements_alert.setWindowTitle("Bad Input")
            N_elements_alert.setText("Entry is not a number, or <= 0")
            N_elements_alert.setInformativeText("Retry")
            N_elements_alert.setIcon(QMessageBox.Icon.Critical)
            N_elements_alert.exec()
        pass

    def setAxis(self, button):
        self.axis = button.text()

    def setNormalized(self, check):
        self.normalize = check

    def ampExcitationToggle(self, check):  # changes visibility of excitation input text edit
        self.ampExcitationList.setRowVisible(self.AmplitudeExcitationList1, check)
        self.ampExcitationList.setRowVisible(self.AmplitudeExcitationList2, check)
        self.ampExcitationList.setRowVisible(self.UniformAmpLabel, not check)

    def updatePhaseValLabel1(self):  # updates the label that informs user of the value of phase shift slider
        real_val1 = self.PhaseSlider1.value()
        self.phaseshift1 = real_val1
        val1 = real_val1 / 100
        if val1:
            self.PhaseValLabel1.setText(f"β_1 = {val1}π")
            self.PhaseValLabel1.setText(f"β_1 = {val1}π")
        else:
            self.PhaseValLabel1.setText(f"β_1 = {val1}")
            self.PhaseValLabel1.setText(f"β_1 = {val1}")

    def updatePhaseValLabel2(self):  # updates the label that informs user of the value of phase shift slider
        real_val2 = self.PhaseSlider2.value()
        self.phaseshift2 = real_val2
        val2 = real_val2 / 100
        if val2:
            self.PhaseValLabel2.setText(f"β_1 = {val2}π")
            self.PhaseValLabel2.setText(f"β_1 = {val2}π")
        else:
            self.PhaseValLabel2.setText(f"β_1 = {val2}")
            self.PhaseValLabel2.setText(f"β_1 = {val2}")

    def resetSlider1(self):
        self.PhaseSlider1.setValue(0)

    def resetSlider2(self):
        self.PhaseSlider2.setValue(0)

    def sweepPhase1(self):
        # HAVE THIS FUNCTION JUST SET THE SLIDER TO LOOP THROUGH SLIDER VALUES THROUGH THREADS
        pass

    def sweepPhase2(self):
        # HAVE THIS FUNCTION JUST SET THE SLIDER TO LOOP THROUGH SLIDER VALUES THROUGH THREADS
        pass

    def updateGraphs(self):
        pass
        #self.TwoDimensional_Plt.graphAF(self.axis, self.N_elements, self.distance, self.phaseshift, self.normalize)
        #self.ThreeDimensional_Plt.graph3DAF(self.axis, self.N_elements, self.distance, self.phaseshift, self.normalize)


class PlanarElevationAzimuthPlot(FigureCanvasQTAgg):
    def __init__(self):
        fig, self.ax = plt.subplots(1, 2, figsize=(10, 10), subplot_kw=dict(polar=True))
        super().__init__(fig)
        fig.tight_layout()

        self.ax[0].grid(True)
        self.ax[1].grid(True)

    def graphAF(self, axis, num_elements, element_distance, phase_shift, normalize):
        self.ax[0].clear()
        self.ax[1].clear()
        # Angular Wave Number * Distance between Elements (k=2π/λ * d=λ/element_distance = 2π/element_distance)
        kd = (2 * np.pi) / element_distance

        # Beta(β) = inter-element phase shift (0 = all elements in phase)
        Beta = np.pi * (phase_shift / 100)

        if axis == "X":
            self.ax[0].set(title='Elevation Plane (XZ-Plane)')
            self.ax[1].set(title='Azimuth Plane (XY-Plane)')
            ElevationPlane_Theta = np.arange(-np.pi, np.pi, 0.01)
            ElevationPlane_Phi = 0
            AzimuthPlane_Theta = np.pi/2
            AzimuthPlane_Phi = np.arange(0, 2 * np.pi, 0.01)
            EP_cos_γ = np.sin(ElevationPlane_Theta) * np.cos(ElevationPlane_Phi)
            AP_cos_γ = np.sin(AzimuthPlane_Theta) * np.cos(AzimuthPlane_Phi)
        elif axis == "Y":
            self.ax[0].set(title='Elevation Plane (YZ-Plane)')
            self.ax[1].set(title='Azimuth Plane (XY-Plane)')
            ElevationPlane_Theta = np.arange(-np.pi, np.pi, 0.01)
            ElevationPlane_Phi = np.pi/2
            AzimuthPlane_Theta = np.pi/2
            AzimuthPlane_Phi = np.arange(0, 2 * np.pi, 0.01)
            EP_cos_γ = np.sin(ElevationPlane_Theta) * np.sin(ElevationPlane_Phi)
            AP_cos_γ = np.sin(AzimuthPlane_Theta) * np.sin(AzimuthPlane_Phi)
        else:
            self.ax[0].set(title='Elevation Plane (Pattern Rotates about Z)')
            self.ax[1].set(title='Azimuth Plane (XY-Plane)')
            ElevationPlane_Theta = np.arange(-np.pi, np.pi, 0.01)
            AzimuthPlane_Theta = np.full(shape=np.arange(0, 2 * np.pi, 0.01).shape, fill_value=np.pi/2)
            EP_cos_γ = np.cos(ElevationPlane_Theta)
            AP_cos_γ = np.cos(AzimuthPlane_Theta)

        # Psi(Ψ) = kd*cos(γ) + β (ELEMENTS ON X AXIS)
        EP_Psi = kd * EP_cos_γ + Beta
        AP_Psi = kd * AP_cos_γ + Beta

        # Array Factor (equal amplitude)
        EP_AF = np.zeros(shape=EP_Psi.shape)  # initializing an array of zeros
        AP_AF = np.zeros(shape=AP_Psi.shape)  # initializing an array of zeros

        for n in range(1, num_elements + 1):
            EP_AF_partial = np.exp(1j * (n - 1) * EP_Psi)  # calculating each element's factor across values of Ψ then concatenating the row
            AP_AF_partial = np.exp(1j * (n - 1) * AP_Psi)  # calculating each element's factor across values of Ψ then concatenating the row
            EP_AF = np.vstack([EP_AF, EP_AF_partial])
            AP_AF = np.vstack([AP_AF, AP_AF_partial])
        EP_AF = abs(np.sum(EP_AF, axis=0))
        AP_AF = abs(np.sum(AP_AF, axis=0))

        if normalize:  # If normalized is checked, normalize the results
            EP_AF = EP_AF/np.max(EP_AF)  # summing matrix columns together and taking absolute value of result, normalizing
            AP_AF = AP_AF/num_elements

        self.ax[0].plot(ElevationPlane_Theta, EP_AF)
        self.ax[1].plot(np.arange(0, 2 * np.pi, 0.01), AP_AF)
        self.draw()


class PlanarFull3DPlot(FigureCanvasQTAgg):
    def __init__(self):
        self.fig, self.ax = plt.subplots(subplot_kw={"projection": "3d"})
        super().__init__(self.fig)
        self.fig.tight_layout()

    def graph3DAF(self, axis, num_elements, element_distance, phase_shift, normalize):
        self.ax.clear()

        Theta = np.arange(0, np.pi, 0.02)  # 0 < θ < 90
        Phi = np.arange(0, 2 * np.pi, 0.02)  # 0 < θ < 360

        # Creating 2D array of Theta and Phi
        THETA, PHI = np.meshgrid(Theta, Phi)

        # Angular Wave Number * Distance between Elements (k=2π/λ * d=λ/element_distance = 2π/element_distance)
        kd = (2 * np.pi) / element_distance

        # Beta(β) = inter-element phase shift (0 = all elements in phase)
        Beta = np.pi * (phase_shift/100)

        if axis == "X":
            cos_γ = np.sin(THETA) * np.cos(PHI)
        elif axis == "Y":
            cos_γ = np.sin(THETA) * np.sin(PHI)
        else:
            cos_γ = np.cos(THETA)

        # Psi(Ψ) = kd*cos(γ) + β
        Psi = kd * cos_γ + Beta

        # Array Factor (equal amplitude)
        AF = np.zeros(shape=Psi.shape)  # initializing an array of zeros

        for n in range(1, num_elements + 1):
            AF_partial = np.exp(1j * (n - 1) * Psi)
            AF = np.add(AF, AF_partial)
        AF = np.abs(AF)

        if normalize:  # If normalized is checked, normalize the results
            AF = AF / np.max(AF)  # summing matrix columns together and taking absolute value of result, normalizing
            self.ax.set_xlim(-1.1, 1.1)
            self.ax.set_ylim(-1.1, 1.1)
            self.ax.set_zlim(-1.1, 1.1)
            #if axis == "Z":
                #self.ax.set_zlim(-1.1, 1.1)
            #else:
                #self.ax.set_zlim(0, 1.1)
        else:
            self.ax.set_xlim(-np.max(AF)-0.1, np.max(AF)+0.1)
            self.ax.set_ylim(-np.max(AF)-0.1, np.max(AF)+0.1)
            #if axis == "Z":
                #pass
            #else:
                #self.ax.set_zlim(0, np.max(AF)+0.1)

        # 3D Plots are in Rectangular coordinates
        X = AF * np.sin(THETA) * np.cos(PHI)
        Y = AF * np.sin(THETA) * np.sin(PHI)
        Z = AF * np.cos(THETA)
        #if axis != "Z":
            #Z[Z < 0] = None  # removes all negative values of Z

        self.ax.plot_surface(X, Y, Z, cmap=matplotlib.colormaps['jet'])

        self.ax.set_xlabel("X")
        self.ax.set_ylabel("Y")
        self.ax.set_zlabel("Z")

        axis = np.linspace(0, 1.1, 2)
        zaxis = np.linspace(0, 1.1, 2)
        axis_zeros = np.zeros(shape=(2,))
        self.ax.set(title='3D Plot')
        self.ax.plot(axis, axis_zeros, zs=0, zdir='z', label='xaxis', color='r')
        self.ax.plot(axis_zeros, axis, zs=0, zdir='z', label='yaxis', color='g')
        self.ax.plot(axis_zeros, axis_zeros, zaxis, label='zaxis', color='b')
        self.ax.legend()
        self.draw()
