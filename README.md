# Mat-Op-Ex Toolkit

Matrix Operator Extraction Toolkit - Tools to obtain the matrix form of common MATLAB linear operations (e.g. trapz, interp) performed over a  1-D or 2-D grid as generated in ndgrid. This solution is based heavily on a wonderful [bit of code](https://www.mathworks.com/matlabcentral/answers/573703-get-interpolation-transfer-relation-matrix-instead-of-interpolated-values) presented by [Bruno Luong](https://www.mathworks.com/matlabcentral/profile/authors/390839), and expanded here to be more generally applied. 

## Description
When performing repeated integration or interpolation of large data-sets over fixed grids, it is often advantageous to access the bare-metal matrix multiplications behind linear functions like interp1, interp2, and trapz. Mat-Op-Ex Toolkit seeks to provide these matrices efficiently, and in a way accessible to the user to build on. 

Mat-Op-Ex functions are designed to be analogous to their conventional counterparts, except that they receive no function data and output a matrix instead. To avoid dealing with high-order tensor operations, physical dimensions are always vectorized within these functions, such that the operator matrix is limited to order N=2. Further details are available in the source code files, and the user should make use of the test scripts available. 

### Function List

  - **Interpolation**
    - `interp1_matrix(X,Xq)` - 1D interp on 1D mesh
    - `interp2_matrix(X,Y,Xq,Yq)` - 1-2D interp on 2D mesh
  - **Integration**
    - `trapz_matrix1(X)` - 1D int. on 1D mesh
    - `trapz_matrix2(X,Y,DIM)` - 1D int. on a 2D mesh 
  - **Auxiliary**
    - `mesh_type` - determines if X/Y inputs are invalid, vectors, *meshgrid* matrices, or *ndgrid* matrices

### Example Test Figures
**(A)** An example spy of the *trapz_matrix2* outout for a 2D grid, applied to either X or Y dimension. 

> > <img src="https://github.com/lynch4815/mat-op-ex/blob/main/figures/test_trapz2_base_spy.png" alt="p1" width="600"/>

**(B)** An example of the stencils used in bi-cubic interpolation in *interp2_matrix*, and the associated deviations from MATLAB interp2 with the cubic method.

> <img src="https://github.com/lynch4815/mat-op-ex/blob/main/figures/test_interp2_base.png" alt="p1" width="600"/>

## Current Status
Mat-Op-Ex is currently tested to support 1D and 2D interpolation (nearest, linear, and cubic) and integration (trapz, simpsons rule) on 1-D and 2-D based grids that are evenly spaced along each axis (i.e. generated with [ndgrid](https://www.mathworks.com/help/matlab/ref/ndgrid.html)). All functions test within machine precision of MATLAB functions, except for the cubic interpolation on the edges of the domain, where a difference in boundary conditions incur a typical error of <1%. 

Because the matrix operators depend heavily on the grid makeup, Mat-Op-Ex functions are less generalized than conventional functions, and expansion to higher dimensions (such as with [interpn](https://www.mathworks.com/help/matlab/ref/interpn.html), are not straightforward. Mat-Op-Ex is similar to [FUNC2MAT](https://www.mathworks.com/matlabcentral/fileexchange/44669-func2mat-convert-linear-function-to-matrix), but far more efficient in generating the matricies. 

## Motivation
This project arose from the need to calculate the stability matrix for a certain set of discretized ODE's where interpolation and spatial integration played a critical role. In addition to furnishing the full matrix operator, the pure matrix form enabled highly-efficient solutions, especially using [gpuArrays](https://www.mathworks.com/help/parallel-computing/run-matlab-functions-on-a-gpu.html) (which to present knowledge, must recompute stencils for interpolation and integration routines, unlike the cpu-based [griddedInterpolant](http://health.ahs.upei.ca/KubiosHRV/MCR/toolbox/matlab/demos/html/griddedInterpolantDemo.html)). Any ideas to improve or expanded on this technique would be greatly welcomed. 

