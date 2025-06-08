# Practice with the deal.II tutorials

Step-23 is a practice for the wave equation.

- In the build case, the mesh is created on an elliptical domain using rectangular minimal units.
- In the merge case, the geometry is created with FreeCAD, and the initial mesh is generated in Gmsh using triangular minimal units.

<br>  
<br>  
<br>  

- In FreeCAD:
  - Use the **Part Design** workbench.
  - Create a **Body** and then a **Sketch**.
  - Draw a reference shape about 5mm × 5mm.
  - Switch to the **Surface Workbench** and use the **Filling** tool.
  - Make sure to create the surface bounded by at least two lines.
  - Export the file via `File → Export` as an **IGES** file.

- In Gmsh software:
  - Open the IGES file by selecting `Open` and navigating to the directory.
  - In Gmsh, go to the **Mesh** menu and select **2D**.
  - Go to `Tools → Options` and check **2D elements faces**.
  - Export the mesh by selecting `File → Export` and save as a **.msh** file.
  - Choose **version 2 ASCII** format when exporting.

<br>  
<br>  
<br>  

## For reference, 

### the official deal.II tutorial:  

https://dealii.org/current/doxygen/deal.II/step_1.html

https://dealii.org/current/doxygen/deal.II/step_2.html

https://dealii.org/current/doxygen/deal.II/step_3.html

https://dealii.org/current/doxygen/deal.II/step_4.html

https://dealii.org/current/doxygen/deal.II/step_23.html

https://dealii.org/current/doxygen/deal.II/step_62.html

### the official deal.II video lectures:

https://www.math.colostate.edu/~bangerth/videos.html

