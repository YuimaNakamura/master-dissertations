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

# 2D Elastic Wave Equation Simulation (deal.II)

This repository implements a 2D elastic wave equation solver using the finite element library deal.II.
The time integration is performed explicitly by the central difference method, and the displacement
field is output as VTU files at each time step for visualization.

---

## Overview

- Uses $Q_1$ (linear) finite elements for the displacement vector $(u_x, u_y)$.
- The spatial domain is a square subdivided 6 times uniformly.
- Two materials are defined by `material_id`, each having different physical properties:
  density $\rho$, and Lamé parameters $\lambda$ and $\mu$.
- Time integration uses an explicit central difference scheme.
- Results are saved in VTU format every 10 time steps for post-processing.

---

## Mathematical Model

### Elastic Wave Equation

The 2D elastic wave equation is given by the PDE:

$$
\rho \frac{\partial^2 \mathbf{u}}{\partial t^2} = \nabla \cdot \boldsymbol{\sigma}
$$

where

- $\mathbf{u} = (u_x, u_y)$ is the displacement vector,
- $\rho$ is the density,
- $\boldsymbol{\sigma}$ is the stress tensor, defined by Hooke's law for linear elasticity with Lamé parameters $\lambda, \mu$.

---

## Numerical Discretization and Time Integration

### Spatial Discretization

Using the finite element method, we approximate the displacement vector $\mathbf{u}$
with basis functions and assemble the stiffness matrix $K$ and mass matrix $M$.

---

### Time Integration: Explicit Central Difference Method

The time stepping scheme for the displacement vector is:

$$
\mathbf{u}^{n+1} = \Delta t^2 M^{-1} (-K \mathbf{u}^n) + 2 \mathbf{u}^n - \mathbf{u}^{n-1}
$$

where

- $\mathbf{u}^n$ is the displacement vector at time step $n$,
- $\Delta t$ is the time step size,
- $M$ is the mass matrix, and $K$ is the stiffness matrix.

---

#### Derivation of the Scheme (Outline)

1. Starting from the semi-discretized equation in space:

$$
M \frac{d^2 \mathbf{u}}{d t^2} + K \mathbf{u} = 0
$$

2. Approximate the second time derivative by the central difference:

$$
\frac{d^2 \mathbf{u}}{d t^2} \approx \frac{\mathbf{u}^{n+1} - 2 \mathbf{u}^n + \mathbf{u}^{n-1}}{\Delta t^2}
$$

3. Substitute into the equation and rearrange:

$$
M \frac{\mathbf{u}^{n+1} - 2 \mathbf{u}^n + \mathbf{u}^{n-1}}{\Delta t^2} + K \mathbf{u}^n = 0
$$

4. Solve for $\mathbf{u}^{n+1}$:

$$
\mathbf{u}^{n+1} = \Delta t^2 M^{-1} (-K \mathbf{u}^n) + 2 \mathbf{u}^n - \mathbf{u}^{n-1}
$$

---

## Code Structure

- The main functionality is implemented in the `ElasticWave` class.
- `setup_system()` initializes finite element spaces and boundary conditions.
- `assemble_system()` assembles the stiffness matrix $K$.
- `assemble_mass_matrix()` assembles the mass matrix $M$ and computes its inverse diagonal entries.
- `initialize_solution()` sets the initial displacement (e.g., Gaussian pulse).
- `time_step()` updates the solution using the explicit central difference scheme.
- `output_results()` writes the solution to VTU files.


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

