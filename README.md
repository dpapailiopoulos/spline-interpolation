# Spline Interpolation

This project implements spline interpolation algorithms for creating smooth curves from a set of 2D points.  
Currently supported methods:  

- **Cubic spline interpolation**  
- **B-spline interpolation**  

The project includes a C++ implementation for performance and a Python script for visualization.  

## Requirements

### C++
- [Eigen3](https://eigen.tuxfamily.org/) (for linear algebra)  
- CMake (for building)  
- A C++17-compatible compiler  

### Python
- [Matplotlib](https://matplotlib.org/) (for visualization)  

## Installation & Build

Clone the repository:
```bash
git clone https://github.com/dpapailiopoulos/spline-interpolation.git
cd spline-interpolation
```
The provided script will take care of building automatically, but manually building is also possible:
```bash
mkdir build && cd build
cmake ..
make
```

## Execution

To run the project, use the ```run.sh``` script:
```bash
./run.sh cubic
```
or 
```bash
./run.sh bspline
```

This will:
1. Build the project (if not already built).
2. Run the selected spline interpolation algorithm.
3. Save the interpolated points to a text file.
4. Visualize the result using Python.