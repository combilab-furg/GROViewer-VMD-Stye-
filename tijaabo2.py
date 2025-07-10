import sys
import os
os.environ["QT_QUICK_BACKEND"] = "software" 
import numpy as np
from PyQt5.QtWidgets import QApplication, QMainWindow, QFileDialog
import pyqtgraph.opengl as gl
import pyqtgraph as pg
import MDAnalysis as mda

class OpenGLGroViewer(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("GRO Molecular Viewer (OpenGL VMD-style)")
        self.setGeometry(100, 100, 1200, 900)
        
        self.view = gl.GLViewWidget()
        self.setCentralWidget(self.view)
        self.view.setCameraPosition(distance=150)
        self.view.opts['center'] = pg.Vector(0, 0, 0)
        
        self.labels = []
        self.bond_threshold = 1.8  # Å
        
        self.load_and_render_gro()

    def get_colors_by_atom_type(self, names):
        """Assign colors based on atom element types (extracted from atom names)."""
        color_map = {
            # Standard CPK coloring scheme (VMD-style)
            'O': (1.0, 0.0, 0.0, 1),    # Red for Oxygen
            'C': (0.2, 0.2, 0.2, 1),    # Dark Gray for Carbon
            'H': (1.0, 1.0, 1.0, 1),    # White for Hydrogen
            'N': (0.0, 0.0, 1.0, 1),    # Blue for Nitrogen
            'P': (1.0, 0.5, 0.0, 1),    # Orange for Phosphorus
            'S': (1.0, 1.0, 0.0, 1),    # Yellow for Sulfur
        }
        # Extract the first character (assumed to be the element) from atom names
        elements = [name[0] if isinstance(name, str) else 'X' for name in names]
        return [color_map.get(elem, (0.8, 0.8, 0.8, 1)) for elem in elements]  # Gray for unknown

    def load_and_render_gro(self):
        try:
            file_path, _ = QFileDialog.getOpenFileName(
                self, "Select .gro or .pdb File", "", "Structure Files (*.gro *.pdb)"
            )
            if not file_path:
                print("No file selected.")
                return
            
            print(f"Loading file: {file_path}")
            u = mda.Universe(file_path)
            print(f"Atoms loaded: {len(u.atoms)}")
            
            positions = u.atoms.positions - u.atoms.positions.mean(axis=0)
            names = u.atoms.names  # Atom names (e.g., 'OW', 'C1', 'H2')
            colors = self.get_colors_by_atom_type(names)  # Color by atom type

            # Render atoms
            for i, pos in enumerate(positions):
                sphere = gl.MeshData.sphere(rows=10, cols=20, radius=0.8)
                atom = gl.GLMeshItem(meshdata=sphere, smooth=True, color=colors[i], shader="shaded")
                atom.translate(*pos)
                self.view.addItem(atom)

            # Render bonds (disable for large systems)
            if len(u.atoms) < 1000:  # Only for small systems
                self.draw_bonds(positions)
            else:
                print("Skipping bonds for large system (too slow).")

            print("Rendering complete.")
            self.showNormal()  # Force window to appear

        except Exception as e:
            print("ERROR:", e)

    def draw_bonds(self, positions):
        """Draw bonds between atoms closer than bond_threshold."""
        n = len(positions)
        for i in range(n):
            for j in range(i+1, n):
                dist = np.linalg.norm(positions[i] - positions[j])
                if dist < self.bond_threshold:
                    pts = np.array([positions[i], positions[j]])
                    bond = gl.GLLinePlotItem(pos=pts, color=(0.7, 0.7, 0.7, 1), width=1.2)
                    self.view.addItem(bond)

if __name__ == "__main__":
    app = QApplication(sys.argv)
    viewer = OpenGLGroViewer()
    viewer.show()
    sys.exit(app.exec_())
