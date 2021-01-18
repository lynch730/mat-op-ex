# Mat-Op-Ex Toolkit

Matrix Operator Extraction Toolkit - Tools to obtain the matrix form of common MATLAB linear operations (e.g. trapz, interp) performed over a  1-D or 2-D grid as generated in ndgrid. This solution is based heavily on a wonderful [bit of code](https://www.mathworks.com/matlabcentral/answers/573703-get-interpolation-transfer-relation-matrix-instead-of-interpolated-values) presented by [Bruno Luong](https://www.mathworks.com/matlabcentral/profile/authors/390839), and expanded here to be more generally applied. 

## Description
When performing repeated integration or interpolation of large datasets over fixed grids, it is often advantegous to access the bare-metal matrix multiplications behind linear functions like interp1, interp2, and trapz. Mat-Op-Ex Toolkit seeks to provide these matricies efficiently, and in a way accessible to the user to build on. 

Mat-Op-Ex functions are designed to be analgous to their conventional counterparts, except that they recieve no function data and output a matrix instead. To avoid dealing with high-order tensor operations, physical dimensions are always vectorized within these functions, such that the operator matrix is limited to order N=2. Further details are available in the source code files, and the user should make use of the test scripts available. 

<a href="url"><img src="https://github.com/lynch4815/mat-op-ex/blob/main/figures/test_interp2_base.png" align="left" height="48" width="48" ></a>

## Current Status
Mat-Op-Ex is currently tested to support 1D and 2D interpolation (nearest, linear, and cubic) and integration (trapz, simpsons rule) on 1-D and 2-D based grids that are evenly spaced along each axis (i.e. generated with [ndgrid](https://www.mathworks.com/help/matlab/ref/ndgrid.html)). All functions test within machine precision of MATLAB functions, except for the cubic interpolation on the edges of the domain, where a difference in boundary conditions incur a typical error of <1%. 

Because the matrix operators depend heavily on the grid makeup, Mat-Op-Ex functions are less generalized than conventional functions, and expansion to higher dimensions (such as with [interpn](https://www.mathworks.com/help/matlab/ref/interpn.html), are not straightforward. Mat-Op-Ex is similar to [FUNC2MAT](https://www.mathworks.com/matlabcentral/fileexchange/44669-func2mat-convert-linear-function-to-matrix), but far more efficient in generating the matricies. 

## Motivation
This project arose from the need to calculate the stability matrix for a certain set of discretized ODE's where interpolation and spatial integration played a critical role. In addition to furnishing the full matrix operator, the pure matrix form enabled highly-efficient solutions, especially using [gpuArrays](https://www.mathworks.com/help/parallel-computing/run-matlab-functions-on-a-gpu.html) (which to present knowledge, must recompute stencils for interpolation and integration routines, unlike the cpu-based [griddedInterpolant](http://health.ahs.upei.ca/KubiosHRV/MCR/toolbox/matlab/demos/html/griddedInterpolantDemo.html)). Any ideas to improve or expanded on this technique would be greatly welcomed. 

