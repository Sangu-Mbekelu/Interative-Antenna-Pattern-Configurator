from PySide6.QtWidgets import QApplication, QMainWindow, QTabWidget
from PySide6.QtGui import QIcon
import sys

from LinearTab import LinearArrayTab
from PlanarTab import PlanarArrayTab
"""
Main Window
Tabs:
- Uniformly-Spaced Linear Element Array Factor Graphing
- Uniformly-Spaced Planar Element Array Factor Graphing
"""
class AFMainWindow(QMainWindow):
    def __init__(self, app):    # Main Window Constructor
        super().__init__()
        self.app = app
        self.setWindowTitle("Array Factor Tool")  # Set Window Title
        self.setWindowIcon(QIcon("Resources/AntennaIcon_Placeholder.png"))  # Set Window Icon

        # Linear Array Widget Class
        linear_pg = LinearArrayTab()

        # Planar Array Widget Class
        planar_pg = PlanarArrayTab()

        # Tabs for Different Linear/Planar Arrays ===============================================================
        tab_menu = QTabWidget()
        tab_menu.addTab(linear_pg, "Linear Phased Array")
        tab_menu.addTab(planar_pg, "Planar Phased Array")

        self.setCentralWidget(tab_menu)


# Defines python Array Pattern Application
AP_App = QApplication(sys.argv)

# Defines Main Widow Interface for Array Pattern Application
Start_Window = AFMainWindow(AP_App)
Start_Window.show()

# Starts event loop - also a blocking function
AP_App.exec()
