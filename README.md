# Practice with the deal.II tutorials

## 🟦 **step-23**

- **Domain**: Unit square or similar simple shape.
- **Mesh**: Structured mesh generated programmatically in deal.II.
- **Elements**: Typically rectangular or quadrilateral.
- **Output**: Solution of the wave equation is written to `.vtu` files.

---

## 🟨 **ellipse**

- **Domain**: Elliptical shape.
- **Mesh**: Created directly in deal.II using a mesh that approximates an ellipse with rectangular/quadrilateral elements.
- **Purpose**: To simulate wave propagation in a smoothly curved domain.
- **Output**: VTU files generated to visualize the wave equation solution over time.

---

## 🟥 **merge**

- **Domain**: Complex geometry created using FreeCAD.
- **Process**:
  1. Design in **FreeCAD**.
  2. Export to `.iges`.
  3. Mesh in **Gmsh** into `.msh` format using triangular elements.
  4. Import into **deal.II** using `GridIn`.
- **Mesh**: Unstructured triangular mesh.
- **Purpose**: To test custom/realistic geometry workflows.
- **Output**: Time-evolved wave equation solution written in `.vtu` files.

---


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

