# Mesh Generation Workflow and Notes (FreeCAD → Gmsh → deal.II)

## Overview

This project involves creating geometries in FreeCAD, generating meshes in Gmsh, and finally importing the `.msh` file into deal.II using the `GridIn` functionality.

## Initial Problem and Solution

### Initial Approach and Issue

Initially, we followed this workflow using the file `simulation_space.FCStd`:

1. Create sketches in **FreeCAD's PartDesign Workbench**.
2. Generate surfaces using the **Surface Workbench**.
3. Select each surface and export it as an **IGES** file.

However, this method resulted in **shared edges being treated as separate entities**.  
For example, even if `Line1` and `Line2` represent the same physical edge, they are assigned different IDs.

When such geometry is processed through Gmsh and imported into deal.II,  
**edges that should be shared are treated as separate**, and deal.II misinterprets them as **external boundaries**.

### Improved Solution

To fix this, we changed the workflow as follows (see `simulation_space_solid.FCStd`):

1. In FreeCAD's **Part Workbench**, combine all surfaces using the `Compound` function.
2. Use the `Solid` creation tool to generate a solid body.
3. Export the geometry as a **STEP file**, ensuring the unit is set to millimeters.
4. Open the STEP file in **Gmsh**.

This method eliminates duplicate edge definitions and ensures **shared edges are correctly recognized**.

## Gmsh Workflow

1. **Open the STEP file** in Gmsh (ensure the unit is millimeters).
2. Use `Geometry → Physical Groups → Add` to assign:
   - **Curves**: for boundary definitions.
   - **Surfaces**: assign `material_id` to each layer.

This process automatically creates a corresponding `.geo` file.

3. Go to `Tools → Options → Geometry → Visibility` to **verify that edge name duplication is resolved**.

4. Under `Tools → Options → Mesh`, make the following changes:
   - Set `2D Algorithm` to `MeshAdapt`.
   - **Uncheck** `Recombine all triangular meshes`.
   - Change `Element size factor` from `1` to `0.1` (for finer mesh resolution).

5. Export the mesh as `simulation1.msh`.

## Importing into deal.II

The generated `simulation1.msh` file can now be loaded into deal.II via the `GridIn` class.  
Each `material_id` assigned via Gmsh's Physical Surfaces is preserved and can be used for assigning material properties or boundary conditions.

---

## Notes

- This setup assumes three layers: `ice`, `saturated_shale_and_silt`, and `gneiss`.
- The `material_id` in deal.II corresponds to the IDs defined in Gmsh's Physical Groups (Surfaces).

<dr>
<dr>
<dr>

📊 Overview of simulation_rootfft_graph_editstep()
This function performs frequency-domain analysis of x-direction displacement data for a specific degree of freedom (DOF) recorded in a ROOT file. It extracts the signal, applies a Fast Fourier Transform (FFT), and visualizes the frequency spectrum.

✅ Main Processing Steps
1. 📂Load ROOT File
- Opens the displacement_data.root file.

- Retrieves the displacement TTree containing simulation data.

2. 🧮Set Up Data Structures
- Prepares variables to receive values from the TTree branches:

- displacement_x: x-direction displacement

- displacement_y: y-direction displacement

- timestep: simulation time step index

- dof_index: degree of freedom index

3. 🔗Connect Branches
- Links the variables to the corresponding TTree branches for data extraction.

4. 🎯Extract Data for a Specific DOF
- Loops through all TTree entries.

- Filters and collects displacement_x values only for dof_index == 1400 into a signal vector.

5. ⚙️Prepare for FFT and Perform Transform
- Initializes arrays for the real and imaginary parts of the signal.

- Uses ROOT’s TVirtualFFT to perform a real-to-complex (R2C) FFT.

6. 📏Define Frequency Axis
- Uses the simulation time step size (delta_t = 1e-4 s) to calculate:

- Sampling rate (e.g., 10,000 Hz)

- Nyquist frequency (half of the sampling rate)

- Sets up a frequency range for plotting the spectrum.

7. 📊Populate FFT Histogram
- Computes the magnitude (amplitude) of each FFT bin from real and imaginary parts.

- Fills a histogram where:

- X-axis = frequency in Hz

- Y-axis = amplitude of displacement

8. 📊Plot and Save Results
- Displays the histogram on a ROOT canvas.

- Saves the plot as a PDF file named fft_simulation_result.pdf.

9. 🧹Cleanup
- Frees dynamically allocated memory used for the FFT.

🎯 Purpose of the Function
- To analyze the frequency content of the x-displacement at a specific DOF (dof_index = 1400) from a time-domain simulation.

- Enables interpretation of wave behavior such as source frequency, reflections, scattering, and boundary effects (e.g., PML) by examining spectral characteristics.

🔧 Requirements
- The ROOT framework must be installed and properly configured.

- The file displacement_data.root must exist and contain a valid displacement TTree.

- The simulation time step size (delta_t) must be known for accurate frequency scaling.