import sys
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
import numpy as np
from genfuncs import *
import SimulatedAnnealingPlanes
import HemisphereMaker
from matplotlib.patches import Patch
from mpl_toolkits.mplot3d import Axes3D
from PyQt5 import QtWidgets, QtCore
import subprocess

def scatter_sphere(ax, x, y, z, r=1, color='b', resolution=15, alpha=0.6):
        """Scatter spheres at x, y, z with radii r (in data units) on 3D axes."""
        def plot_sphere(center, radius):
            u = np.linspace(0, 2 * np.pi, resolution)
            v = np.linspace(0, np.pi, resolution)
            x_s = center[0] + radius * np.outer(np.cos(u), np.sin(v))
            y_s = center[1] + radius * np.outer(np.sin(u), np.sin(v))
            z_s = center[2] + radius * np.outer(np.ones_like(u), np.cos(v))
            return ax.plot_surface(x_s, y_s, z_s, color=color, alpha=alpha, linewidth=0, shade=True)

        # Handle broadcasting of inputs
        x = np.atleast_1d(x)
        y = np.atleast_1d(y)
        z = np.atleast_1d(z)
        r = np.broadcast_to(r, x.shape)

        for xi, yi, zi, ri in zip(x, y, z, r):
            plot_sphere((xi, yi, zi), ri)
setattr(Axes3D, 'scatter_sphere', scatter_sphere)
# --------------------------------------------------
# Matplotlib 3D Canvas for displaying lattice
# --------------------------------------------------
class LatticeCanvas(FigureCanvas):
    def __init__(self, parent=None):
        self.fig = Figure()
        self.ax = self.fig.add_subplot(111, projection="3d")
        self.ax.set_axis_off()  # Hide axes
        super().__init__(self.fig)

    def plot_latt(self, lattice, nuc_scale = 1, show_bonds=False):
        self.ax.clear()
        self.ax.set_axis_off()
        l = {}
        for atom in lattice.atoms:
            # Draw atom as scatter point with color + size
            x,y,z = lattice.get_atom_cart_coords(atom)
            #change res to help improve performance. I know its kinda trash, but if you want a fast cif visualizer, use vista.
            res = max(8,250/len(lattice.atoms)) #set min
            res = min(res,50) #set max
            res = int(res) #ensure int

            self.ax.scatter_sphere(
                x,y,z,
                r=atom.radius * nuc_scale,
                color=ELEMENT_COLORS[atom.element],
                resolution=res
            )
            if not(atom.element in l):
                l[atom.element] = ELEMENT_COLORS[atom.element]
            self.ax.set_aspect('equal', 'box')
        if show_bonds:
            for bond in lattice.bonds:
                x1, y1, z1 = lattice.get_atom_cart_coords(lattice.atoms[bond[0]])
                x2, y2, z2 = lattice.get_atom_cart_coords(lattice.atoms[bond[1]])
                self.ax.plot((x1,x2), (y1,y2), (z1,z2), color="k", linewidth=3, alpha=0.9)

        handles = [Patch(color=color, label=label) for label, color in l.items()]
        self.ax.legend(handles=handles, title='Atoms', loc='upper right')
        self.draw()

# --------------------------
# -- Widgets
# --------------------------
class FloatSlider(QtWidgets.QWidget):
    """Composite widget: a slider + line edit to control a float value."""

    value_changed = QtCore.pyqtSignal(float)

    def __init__(self, label, minimum=0.0, maximum=1.0, step=0.01, default=None, parent=None):
        super().__init__(parent)
        self.min, self.max, self.step = minimum, maximum, step
        self.scale = int(round(1.0 / step))  # scale factor to map float to int slider positions
        if default is None:
            default = (minimum + maximum) / 2.0

        layout = QtWidgets.QHBoxLayout(self)

        # Label on the left
        self.name_label = QtWidgets.QLabel(label)
        self.name_label.setFixedWidth(140)
        layout.addWidget(self.name_label)

        # Slider in the middle
        self.slider = QtWidgets.QSlider(QtCore.Qt.Horizontal)
        self.slider.setMinimum(int(round(minimum * self.scale)))
        self.slider.setMaximum(int(round(maximum * self.scale)))
        self.slider.setValue(int(round(default * self.scale)))
        layout.addWidget(self.slider, stretch=1)

        # Text box for precise entry
        self.value_display = QtWidgets.QLineEdit(f"{default:.6g}")
        self.value_display.setFixedWidth(90)
        self.value_display.setAlignment(QtCore.Qt.AlignRight)
        layout.addWidget(self.value_display)

        # Connect signals
        self.slider.valueChanged.connect(self._on_slider)
        self.value_display.editingFinished.connect(self._on_edit)

    def _on_slider(self, int_val):
        """Update text box when slider moves."""
        f = int_val / self.scale
        self.value_display.setText(f"{f:.6g}")
        self.value_changed.emit(f)

    def _on_edit(self):
        """Update slider when text box is edited."""
        try:
            f = float(self.value_display.text())
        except ValueError:
            f = self.slider.value() / self.scale
        f = max(self.min, min(self.max, f))  # clamp
        self.slider.setValue(int(round(f * self.scale)))
        self.value_display.setText(f"{f:.6g}")
        self.value_changed.emit(f)

    def get(self):
        """Return the current float value."""
        return self.slider.value() / self.scale

    def set(self, v):
        """Set the slider/textbox value programmatically."""
        v = max(self.min, min(self.max, v))
        self.slider.setValue(int(round(v * self.scale)))
        self.value_display.setText(f"{v:.6g}")


class IntSpinBox(QtWidgets.QWidget):
    """Composite widget: label + QSpinBox for integer parameters."""

    value_changed = QtCore.pyqtSignal(int)

    def __init__(self, label, minimum=1, maximum=1000, default=10, parent=None):
        super().__init__(parent)
        layout = QtWidgets.QHBoxLayout(self)

        # Label
        self.name_label = QtWidgets.QLabel(label)
        self.name_label.setFixedWidth(140)
        layout.addWidget(self.name_label)

        # Spin box
        self.spin = QtWidgets.QSpinBox()
        self.spin.setRange(minimum, maximum)
        self.spin.setValue(default)
        layout.addWidget(self.spin)

        # Connect
        self.spin.valueChanged.connect(self.value_changed.emit)

    def get(self):
        return self.spin.value()

    def set(self, v):
        self.spin.setValue(bool(v))

# --------------------------
# -- Main Window
# --------------------------
class MainWindow(QtWidgets.QMainWindow):
    """The main application window."""

    def __init__(self):
        super().__init__()
        self.setWindowTitle("Crystal Lattice GUI")
        self.data_path = 'CustomData.dat'
        # self.resize(1150, 780)

        # Store unit mode (internal is always mm)
        self.unit_is_inches = False

        # ---------------- Central Layout ----------------
        central = QtWidgets.QWidget()
        self.setCentralWidget(central)
        layout = QtWidgets.QHBoxLayout(central)

        # Left panel: controls
        controls = QtWidgets.QVBoxLayout()
        layout.addLayout(controls, stretch=0)

        
        # -- File input
        in_widget = QtWidgets.QWidget()
        in_layout = QtWidgets.QHBoxLayout(in_widget)
        self.infile_edit = QtWidgets.QLineEdit()
        in_btn = QtWidgets.QPushButton("Browse")
        in_btn.clicked.connect(self._browse_infile)
        in_layout.addWidget(QtWidgets.QLabel("Input file:"))
        in_layout.addWidget(self.infile_edit)
        in_layout.addWidget(in_btn)
        controls.addWidget(in_widget)

        # -- File output
        self.out_widget = QtWidgets.QWidget()
        self.out_layout = QtWidgets.QHBoxLayout(self.out_widget)
        self.outfile_edit = QtWidgets.QLineEdit()
        self.out_btn = QtWidgets.QPushButton("Browse")
        self.out_btn.clicked.connect(self._browse_outfile)
        self.out_layout.addWidget(QtWidgets.QLabel("Output path:"))
        self.out_layout.addWidget(self.outfile_edit)
        self.out_layout.addWidget(self.out_btn)
        controls.addWidget(self.out_widget)
        self.out_widget.hide()


        # -- parameters
        # Format: name: (min, max, step, default, pregenerate, advanced_only)
        self.param_defs = { #note, they are displayed in this order
            "nx": (1, 10.0, 1, 1.0, False, False),
            "ny": (1, 10.0, 1, 1.0, False, False),
            "nz": (1, 10.0, 1, 1.0, False, False),
            "Curve Resolution": (3, 500, 1, 50, False, False),
            "Buildplate Length": (10, 1000, 10, 210, False, False),
            "Buildplate Width": (10, 1000, 10, 250, False, False),
            "Radius": (0.1, 30.0, 0.1, 15.0, False, False),
            "Nucleus Scale": (.01, 1.0, .01, 1.0, False, False),
            "Bond Range": (0.0, 1, 0.01, 0.1, False, True),
            "Tolerance": (0.0, 1.0, 0.01, 0.1, False, True),
            "Cut Depth Ratio": (0.01, 1.0, 0.01, 0.25, False, True),
            "Cut Radius Ratio": (0.01, 1.0, 0.01, 0.16, False, True),
            "Trapezoid Side Ratio": (0.01, 1.0, 0.01, 0.33, False, True),
            "Trapezoid Depth Ratio": (0.01, 1.0, 0.01, 0.2, False, True),
            "Trapezoid Oblong Ratio": (0.01, 1.0, 0.01, 0.9, False, True),
            "Overhang Angle": (0.0, 90.0, 0.5, 45.0, False, True)
        }

        # -- Container for parameters
        self.param_container = QtWidgets.QWidget()
        self.param_layout = QtWidgets.QVBoxLayout(self.param_container)
        controls.addWidget(self.param_container)

        # store slider widgets with metadata
        self.sliders = {}
        self.slider_widgets = []  # (widget,pre_gen, is_advanced)

        for name,(mn,mx,step,default,pre_gen,is_adv) in self.param_defs.items():
            if step%1 != 0:
                s = FloatSlider(name,float(mn),float(mx),float(step),float(default))
                self.sliders[name] = s
                w = QtWidgets.QWidget()
                l = QtWidgets.QVBoxLayout(w)
                l.setContentsMargins(0,0,0,0)
                l.addWidget(s)
                self.slider_widgets.append((w,pre_gen,is_adv))
            else:
                c = IntSpinBox(name,int(mn),int(mx),int(default))
                self.sliders[name] = c
                w = QtWidgets.QWidget()
                l = QtWidgets.QVBoxLayout(w)
                l.setContentsMargins(0,0,0,0)
                l.addWidget(c)
                self.slider_widgets.append((w,pre_gen,is_adv))

        # initially apply beginner layout
        self.generated = False #needed for apply slider layout
        self._apply_slider_layout(advanced=False)


        # -- Bonds checkbox
        self.bonds_checkbox = QtWidgets.QCheckBox("Show bonds")
        self.bonds_checkbox.setVisible(False)
        self.bonds_checkbox.stateChanged.connect(self._toggle_bonds)
        controls.addWidget(self.bonds_checkbox)

        # -- Beginner/Advanced toggle
        self.mode_checkbox = QtWidgets.QCheckBox("Advanced mode")
        self.mode_checkbox.stateChanged.connect(self._on_mode_change)
        self.mode_checkbox.setVisible(False)
        controls.addWidget(self.mode_checkbox)

        # -- Separate Unique atoms and bonds toggle
        self.separate_checkbox = QtWidgets.QCheckBox("Separate atoms types and bonds")
        self.separate_checkbox.stateChanged.connect(self._separate_atoms)
        self.separate_checkbox.setVisible(False)
        controls.addWidget(self.separate_checkbox)

        # -- Units toggle
        self.unit_widget = QtWidgets.QWidget()  # container so we can hide/show as one block
        unit_layout = QtWidgets.QHBoxLayout(self.unit_widget)
        unit_layout.addWidget(QtWidgets.QLabel("Units:"))
        self.mm_radio = QtWidgets.QRadioButton("mm")
        self.in_radio = QtWidgets.QRadioButton("inches")
        self.mm_radio.setChecked(True)
        unit_layout.addWidget(self.mm_radio)
        unit_layout.addWidget(self.in_radio)
        controls.addWidget(self.unit_widget)
        self.unit_widget.hide()  # hidden at startup

        # -- Action buttons
        btns = QtWidgets.QHBoxLayout()
        self.gen_btn = QtWidgets.QPushButton("Generate")
        self.exp_btn = QtWidgets.QPushButton("Export")
        self.gen_btn.clicked.connect(self._on_generate)
        self.exp_btn.clicked.connect(self._on_export)
        self.exp_btn.setEnabled(False)
        btns.addWidget(self.gen_btn); btns.addWidget(self.exp_btn)
        controls.addLayout(btns)

        # -- Log box
        self.log_box = QtWidgets.QTextEdit()
        self.log_box.setReadOnly(True)
        # self.log_box.setMinimumHeight(150)
        controls.addWidget(QtWidgets.QLabel("Log:"))
        controls.addWidget(self.log_box)


        # add version label
        controls.addWidget(QtWidgets.QLabel("Version 1.1"))
        controls.addStretch()

        # Right panel: lattice viewer
        self.canvas = LatticeCanvas()
        layout.addWidget(self.canvas, stretch=1)

        # Apply initial mode visibility
        self._on_mode_change()
        

    # --------------------------
    # Helper methods
    # --------------------------
    def log_message(self, text, error=False):
        """Append a message to the log window (red if error)."""
        color = "red" if error else "black"
        self.log_box.append(f"<span style='color:{color}'>{text}</span>")
        self.log_box.ensureCursorVisible()
        QtWidgets.QApplication.processEvents()  # forces repaint & log update

    def _browse_infile(self):
        """Open a file dialog restricted to .cif files."""
        path, _ = QtWidgets.QFileDialog.getOpenFileName(
            self,
            "Select input CIF file",
            filter="CIF files (*.cif)"
        )
        if path:
            self.infile_edit.setText(path)

    def _browse_outfile(self):
        """Open a folder dialog for choosing an export directory."""
        path = QtWidgets.QFileDialog.getExistingDirectory(
            self,
            "Select export folder"
        )
        if path:
            # Store folder path
            self.outfile_edit.setText(path)

    def _on_mode_change(self):
        """Show/hide advanced parameters depending on toggle state."""
        self._apply_slider_layout(self.mode_checkbox.isChecked())

    def _separate_atoms(self):
        """Show/hide advanced parameters depending on toggle state."""
        # self._apply_slider_layout(self.separate_checkbox.isChecked())
        pass

    def _on_unit_change(self):
        """Switch displayed slider values between mm and inches."""
        self.unit_is_inches = self.in_radio.isChecked()

    def _collect_params(self):
        """Gather all parameter values into a dict, converting to mm."""
        params = {}
        for name in self.param_defs:
            val = self.sliders[name].get()
            if (name == 'Radius' or name == 'Tolerance') and self.unit_is_inches:  # convert to mm
                val *= 25.4
            params[name] = val
        return params

    def _on_generate(self):
        infile = self.infile_edit.text().strip()
        if len(infile) == 0:
            self.log_message("Select input cif file.",True)
        elif infile[-4:] != '.cif':
            self.log_message("Ensure input file is cif file.",True)
        else:
            if not self.generated: #make changes for generating STL
                self.generated = True
                self.gen_btn.setText("Regenerate")
                self.exp_btn.setEnabled(True)
                self.mode_checkbox.setVisible(True)
                self.separate_checkbox.setVisible(True)
                self.unit_widget.show()
                self.out_widget.show()
                self.bonds_checkbox.setVisible(True)
                self._apply_slider_layout(advanced=False)
                self.out_btn.setVisible(True)
                self.outfile_edit.setVisible(True)
            
            
            
            params = self._collect_params()
            nx,ny,nz = params['nx'],params['ny'],params['nz']
            self.lattice = Lattice(infile,nx,ny,nz,bond_range=params['Bond Range']) #parce cif file and create a lattice object with all spacial and atomic data (see genfuncs.py)
            self.log_message(f"{self.lattice.lattice_name} lattice generated. {len(self.lattice.atoms)} atoms. {len(self.lattice.bonds)} bonds.",False)
            if not self.lattice.all_nodes_connected():
                self.log_message("Warning: Not all atoms are connected. Consider increasing bond range.",True)
            self.canvas.plot_latt(self.lattice, nuc_scale=params["Nucleus Scale"], show_bonds=self.bonds_checkbox.isChecked())

    def _on_export(self):
        """Handler for Export button."""
        openscad_path = self._find_OpenSCAD()
        params = self._collect_params()
        outfile = self.outfile_edit.text().strip()
        if not outfile:
            self.log_message("No output path selected", error=True)
            return
        cut_planes = []
        for atom in self.lattice.atoms:
            optimal_plane_normal, _ = SimulatedAnnealingPlanes.simulated_anneal_atoms(atom.bonds)
            cut_planes.append(optimal_plane_normal)
        ps = [
            params["Curve Resolution"], #number of fragments for each 3D printed piece. basically 3D printed circle resolution
            params["Radius"],#max radius size in mm
            params["Nucleus Scale"],# scale from 0 to 1 to decrease nucleus size in order to see structure more clearly
            params["Tolerance"],#add tolerance to male parts
            params["Cut Depth Ratio"],# Depth of the cylindrical cut for bonds
            params["Cut Radius Ratio"],# Radius of the cutting cylinder
            params["Trapezoid Side Ratio"],#ratio of half long trap leg to radius
            params["Trapezoid Depth Ratio"],# ratio of trap prism cut depth to radius
            params["Trapezoid Oblong Ratio"],#ratio of height of short trap leg to long trap leg
            params["Overhang Angle"]#degrees, internal pryamid that prevents supports in female hemisphere. Should be >45 and <90
        ]
        HemisphereMaker.makeSTL(outfile,self.lattice.atoms,cut_planes,plate_dim=(params['Buildplate Length'],params['Buildplate Width']),parameters=ps,log=self, sep = self.separate_checkbox.isChecked(),openscad_cmd = openscad_path)

    def _toggle_bonds(self, state):
        """Toggle bond drawing in the 3D lattice."""
        params = self._collect_params()
        self.canvas.plot_latt(self.lattice, nuc_scale=params["Nucleus Scale"], show_bonds=self.bonds_checkbox.isChecked())
    
    def _apply_slider_layout(self, advanced):
        """Rebuild the parameter layout in beginner or advanced mode."""
        #this is pretty janky and probably fairly easy to break...
        # clear old layout
        for i in reversed(range(self.param_layout.count())):
            item = self.param_layout.takeAt(i)
            if item.widget():
                item.widget().setParent(None)
        if not self.generated:#single column
            vbox = QtWidgets.QVBoxLayout()
            for w, pre_gen, is_adv in self.slider_widgets:
                if not(is_adv) and pre_gen:
                    vbox.addWidget(w)
            vbox.addStretch()
            self.param_layout.addLayout(vbox)
        else:#double column
            if not advanced:
                # Beginner mode: 
                grid = QtWidgets.QGridLayout()
                row,col = 0,0
                for w, pre_gen, is_adv in self.slider_widgets:
                    if not is_adv:
                        grid.addWidget(w, row, col)
                        col += 1
                        if col == 2:  # wrap to next row
                            col = 0
                            row += 1
                self.param_layout.addLayout(grid)
            else:
                # Advanced mode: show everything in two columns using a grid
                grid = QtWidgets.QGridLayout()
                row,col = 0,0
                for w, _, is_adv in self.slider_widgets:
                    grid.addWidget(w, row, col)
                    col += 1
                    if col == 2:  # wrap to next row
                        col = 0
                        row += 1
                self.param_layout.addLayout(grid)
    
    def _render_time_estimate(self):
        atom_num = len(self.lattice.atoms)
        params = self._collect_params()
        res = params["Curve Resolution"]
        a = 0.105#roughly fitted parameters
        b = 0.02405#roughly fitted parameters
        c = -0.125#roughly fitted parameters
        return atom_num*(a*np.exp(b*res)+c)

    def _find_OpenSCAD(self): 
        #this file is meant to allow for recursive calling in order to ensure the correct OpenScad.exe is found
        
        openscad_path = self._get_datafile_path(self.data_path) #get path from custom data file

        try:
            result = subprocess.run(
                [openscad_path, "--version"],  # just checking to make sure OpenSCAD exists
                capture_output=True,
                text=True,
                check=True
            )
            print(result.stderr.strip())
            return openscad_path
        except FileNotFoundError:
            #run popup to find OpenSCAD
            QtWidgets.QMessageBox.information(
                self, "OpenSCAD Not Found", f"Naviagate to OpenSCAD location in subsequent popup."
            )
            openscad_path2, _ = QtWidgets.QFileDialog.getOpenFileName(
                self,
                "Select OpenSCAD.exe file",
                "",
                "EXE files (*.exe);;All files (*)"
            )

            if not openscad_path2: #if close prompt without giving path
                return self._find_OpenSCAD()
            

            # Add OpenSCAD file to data file
            with open(self.data_path, "r") as f:
                lines = f.readlines()
            lines[0] = lines[0].strip().split()[0] + ' ' + openscad_path2 + '\n'

            # Write back
            with open(self.data_path, "w") as f:
                f.writelines(lines)

            # Show success message
            QtWidgets.QMessageBox.information(
                self, "Success", f"OpenSCAD location saved successfully"
            )
            return self._find_OpenSCAD()#repeat to ensure the located file is correct

    def _get_datafile_path(self,data_path):
        with open(data_path, "r") as f:
            lines = f.readlines()
        line = lines[0]
        if len(line)>14:
            return line[14:].strip()
        else:
            print(f'Failed to read data file line:{lines[0]}')
            return
        
    def _register_render_time(self,res,t):
        # Add OpenSCAD file to data file
        with open(self.data_path, "r") as f:
            lines = f.readlines()
        lines.append(f'{res},{t},{len(self.lattice.atoms)},{len(self.lattice.bonds)}\n')

        # Write back
        with open(self.data_path, "w") as f:
            f.writelines(lines)

# --------------------------
# -- Run
# --------------------------
app = QtWidgets.QApplication(sys.argv)
win = MainWindow()
win.show()
sys.exit(app.exec_())


# command to make debugging
# python -m PyInstaller LatticetoSTL.py --onefile --windowed -i "C:\Users\blake\Downloads\LogoICO.ico" -d all --console

# command to make reg version
# python -m PyInstaller LatticetoSTL.py --onefile --windowed -i "C:\Users\blake\Downloads\LogoICO.ico"