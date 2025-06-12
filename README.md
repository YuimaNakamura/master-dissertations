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

## 🟦 **elastic**

# 2D Elastic Wave Equation Solver using deal.II (Explicit Central Difference Scheme)

This project implements a 2D elastic wave simulation using the finite element library [deal.II](https://www.dealii.org/), with an explicit central difference scheme for time integration.

## 📌 Overview

We solve the linear **elastic wave equation** in two dimensions using:
- Finite element method (FEM) for **spatial discretization**.
- Central difference method (explicit) for **temporal integration**.

The domain may consist of two adjacent subdomains (e.g., materials A and B) with **different material properties** (density, Lamé constants).

## 📐 Mathematical Model

We consider the second-order hyperbolic system:

\[
\rho(\mathbf{x}) \frac{\partial^2 \mathbf{u}}{\partial t^2} = \nabla \cdot \boldsymbol{\sigma}(\mathbf{u}) + \mathbf{f}(\mathbf{x}, t)
\]

- \(\mathbf{u}(\mathbf{x},t)\): displacement vector  
- \(\rho(\mathbf{x})\): density (possibly spatially varying)  
- \(\mathbf{f}(\mathbf{x},t)\): external body force (e.g., Gaussian pulse)  
- \(\boldsymbol{\sigma}(\mathbf{u})\): stress tensor defined by Hooke’s law

### Constitutive Law (Linear Isotropic Elasticity)

\[
\boldsymbol{\sigma} = \lambda (\nabla \cdot \mathbf{u}) \mathbf{I} + 2\mu \, \boldsymbol{\varepsilon}(\mathbf{u})
\]
\[
\boldsymbol{\varepsilon}(\mathbf{u}) = \frac{1}{2}(\nabla \mathbf{u} + \nabla \mathbf{u}^T)
\]

- \(\lambda, \mu\): Lamé parameters (can vary between subdomains)

## 🧮 Weak Form and Discretization

### Spatial Discretization (FEM)

The variational (weak) form after multiplying by test functions and integrating:

\[
\int_\Omega \rho \frac{\partial^2 \mathbf{u}_h}{\partial t^2} \cdot \mathbf{v}_h \, dx + \int_\Omega \boldsymbol{\sigma}(\mathbf{u}_h) : \nabla \mathbf{v}_h \, dx = \int_\Omega \mathbf{f} \cdot \mathbf{v}_h \, dx
\]

This leads to a semi-discrete matrix system:

\[
\mathbf{M} \ddot{\mathbf{u}}(t) + \mathbf{K} \mathbf{u}(t) = \mathbf{F}(t)
\]

- \(\mathbf{M}\): consistent mass matrix  
- \(\mathbf{K}\): stiffness matrix  
- \(\mathbf{u}(t)\): global displacement vector  
- \(\mathbf{F}(t)\): external force vector

## ⏱️ Time Discretization (Explicit Central Difference Method)

We use a second-order central difference method for time integration:

### Second-order time derivative:

\[
\ddot{\mathbf{u}}^n \approx \frac{\mathbf{u}^{n+1} - 2\mathbf{u}^n + \mathbf{u}^{n-1}}{\Delta t^2}
\]

Substitute into the semi-discrete system:

\[
\mathbf{M} \left( \frac{\mathbf{u}^{n+1} - 2\mathbf{u}^n + \mathbf{u}^{n-1}}{\Delta t^2} \right) + \mathbf{K} \mathbf{u}^n = \mathbf{F}^n
\]

Rearranging gives the explicit update formula:

\[
\boxed{
\mathbf{u}^{n+1} = \Delta t^2 \, \mathbf{M}^{-1} (\mathbf{F}^n - \mathbf{K} \mathbf{u}^n) + 2\mathbf{u}^n - \mathbf{u}^{n-1}
}
\]

This scheme requires only the current and previous displacement vectors, making it efficient for time stepping.

## ⚙️ Code Structure

- `elastic_wave.cc`: Main simulation code using deal.II.
- `CMakeLists.txt`: Build configuration.
- Output files: Time-stepped VTU files for visualization (e.g., `solution_0001.vtu`).

### Key Components
- `ElasticWaveProblem` class:
  - `setup_system()`: Mesh generation, DoF handling, mass/stiffness matrix setup
  - `assemble_system()`: Assembles \(\mathbf{M}\), \(\mathbf{K}\)
  - `solve_timestep()`: Explicit time stepping via above formula
  - `output_results()`: VTU output using DataOut

## ▶️ How to Build and Run

```bash
mkdir build
cd build
cmake ..
make
./elastic_wave


---


<br>  
<br>  
<br>  

- 🔵In FreeCAD:
  - Use the **Part Design** workbench.
  - Create a **Body** and then a **Sketch**.
  - Draw a reference shape about 5mm × 5mm.
  - Switch to the **Surface Workbench** and use the **Filling** tool.
  - Make sure to create the surface bounded by at least two lines.
  - Export the file via `File → Export` as an **IGES** file.

- 🟢In Gmsh software:
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

