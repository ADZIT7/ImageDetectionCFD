# **ImageDetectionCFD**

A specialized Computational Fluid Dynamics (CFD) toolkit that bridges the gap between digital images and aerodynamic simulation. This project allows users to extract geometric shapes from images using Computer Vision and perform flow analysis using the **Source Panel Method** and **Navier-Stokes solvers**.

## **Features**

* **Geometry Extraction:** Uses OpenCV and Canny edge detection to extract $(x, y)$ coordinates from images (e.g., airfoils, cylinders, or custom shapes).  
* **Source Panel Method:** Implements a 2D potential flow solver to calculate:  
  * Pressure distribution ($C\_p$) across the body surface.  
  * Lift and Drag coefficients.  
  * Kutta condition enforcement for lifting surfaces.  
* **Navier-Stokes Solver:** Includes a 2D incompressible flow solver (Lid-Driven Cavity) to model viscous flow patterns.  
* **Visualization:** Generates high-quality streamline plots, pressure contours, and velocity fields using Matplotlib.

## **Project Structure**

* shapesdetection.py: Processes input images (.png, .jpg) to extract edge coordinates.  
* ObjectDetectionAirFoil.py: The main solver. Interpolates shape data and runs the Source Panel Method simulation.  
* CFD.py: Core logic for potential flow calculations and streamline generation.  
* multi.py: A finite-difference solver for the Incompressible Navier-Stokes equations.

## **Installation**

1. **Clone the repository:**  
   ```
   git clone https://github.com/ADZIT7/ImageDetectionCFD.git
   cd ImageDetectionCFD
   ```
   
2. **Create and activate a virtual environment**

   Using a virtual environment is strongly recommended to isolate dependencies.

   **macOS / Linux**
   ```
   python3 -m venv venv
   source venv/bin/activate
   ```
  
   **Windows (PowerShell)**
   ```
   python -m venv venv
   venv\Scripts\activate
   ```
  
   Once activated, your terminal prompt should show `(venv)`.

3. **Install dependencies:**  
   ```
   pip install \-r requirements.txt
   ```

## **Usage**

### **1\. Shape Extraction**

Extract coordinates from an image (default looks for circle.png):

python shapesdetection.py

### **2\. Run CFD Simulation**

Ensure your coordinate file (e.g., square2.dat) is referenced in ObjectDetectionAirFoil.py, then run:

python ObjectDetectionAirFoil.py

## **Mathematical Background**

The potential flow solver represents the body as a series of $N$ panels. It solves a system of linear equations to find source strengths ($\\lambda$) and circulation ($\\Gamma$) such that the flow is tangent to the surface at every control point:

$$ \\sum\_{j=1}^{N} \\lambda\_j A\_{ij} \+ \\vec{V}\_{\\infty} \\cdot \\vec{n}\_i \= 0 $$

The Navier-Stokes solver (multi.py) utilizes a Pressure-Poisson equation to ensure mass conservation in incompressible flow.

## **License**

Distributed under an All Rights Reserved License. See `LICENSE` for more information. Contact me if you’re interested in purchasing a license for commercial use.
