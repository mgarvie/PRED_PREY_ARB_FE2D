## Repository:
 
PRED_PREY_ARB_FE2D is a collection of MATLAB codes using the finite element method for simulating Predator-Prey Interactions in 2D.

## General Description:

PRED_PREY_ARB is a collection of simple MATLAB routines using the finite element method for simulating the dynamics of predator-prey interactions modelled by a nonlinear reaction-diffusion system. Unlike FD2D the systems are solved on domains of arbitrary shape using general boundary conditions. The collection of 2D codes are called FE2D. 

FE2D is a collection of MATLAB routines using the finite element method for simulating the dynamics of predator-prey interactions in two space dimensions and time. The codes generalize my earlier predator-prey codes solved on the square (see FD2D) to the situation of arbitrary shaped domains with general boundary conditions (functions depending on space and time of Neumann, Dirichlet, Robin, Periodic, or mixed type).

The MATLAB code is mostly self explanatory, with the names of variables and parameters corresponding to the symbols used in the finite element method described in the paper referenced below. Copies of the MATLAB codes are freely available via the links below.

The code employs the sparse matrix facilities of MATLAB when solving the linear systems, which provides advantages in both matrix storage and computation time. The code is vectorized to minimize the number of "for-loops" and conditional "if-then-else" statements, which again helps speed up the computations.

The linear systems are solved using MATLAB's built in function gmres.m. The gmres.m algorithm in MATLAB requires a number of input arguments. For some simulations we found it acceptable to use no "restarts" of the iterative method, or, preconditioners, and a tolerance for the relative error of 1E-6. In practise, the user will need to experiment with the restart value, tolerance, and "maxit" (maximum number of outer iterations) in order to achieve satisfactory rates of convergence of the gmres.m function. For definitions of these input arguments (and others) see the description in MATLAB's Help pages. We remark that a pure C or FORTRAN code is likely to be faster than our codes, but with the disadvantage of much greater complexity and length.

The user is prompted for all the necessary parameters, time-step, initial data functions, and boundary data functions. However, the grid must be supplied by the user (further details given below).

The reaction-diffusion system models spatially extended predator-prey interactions of Rosenzweig-MacArthur form, where u(x,t) and v(x,t) are the prey and predator densities at time t and vector position x, and Delta is the usual Laplacian operator, and the parameters alpha, beta, gamma, and delta are strictly positive.

## Grid generation

Before the finite element codes can be run in MATLAB it is necessary to triangulate the region using an appropriate unstructured grid generator in 2D. As a starting point, we assume that the user has provided a list of node coordinates, ordered clock-wise or counter-clockwise, corresponding to the boundary of the computational domain. For simple geometries this is easily worked out by hand, but for more complicated domains the use of a drawing package may be advisable. The grid generator then needs to output two arrays that define the mesh, namely

    t, a list of node indices corresponding to M triangles, given by an 3-by-M array;
    p, a list of N node coordinates, given by an 2-by-N array;

We found it convenient to use an unstructured mesh generator (Mesh2d v24) in MATLAB freely available from http://www.mathworks.com/matlabcentral/fileexchange/25555-mesh2d-automatic-mesh-generation". 

## The initial data

The user is prompted for the initial data, which can be functions depending on space, or just constants. Unlike with the simpler PRED_PREY_SIM codes no special format is required for entering the initial data functions.

An example with an acceptable format is the following:

        >> Enter initial data function u0(x,y)  0.2*exp(-(x^2+y^2))
        >> Enter initial predator function v0(x,y)  1.0
      

We can also define functions in a piecewise fashion. For example, with Omega=[0,100]x[0,100], in order to choose an initial predator density of 0.2 within a circle with radius 5 and center (50,50), and a density of 0 elsewhere on Omega we input the following:
      
        >> Enter initial data function u0(x,y) 0.2*((x-50)^2+(y-50)^2<25)
      
## The boundary conditions

The boundary of the whole domain is denoted Gamma. In the mixed boundary condition cases we assume Gamma is the union of two pieces Gamma1 and Gamma2. The code allows for the following boundary conditions:

    Pure Neumann boundary conditions on Gamma.
    Pure Dirichlet boundary conditions on Gamma.
    Pure Robin boundary conditions on Gamma.
    Mixed boundary conditions with a Dirichlet condition on Gamma1 and a Neumann condition on Gamma2.
    Mixed boundary conditions with a Robin condition on Gamma1 and a Neumann condition on Gamma2.
    Periodic boundary conditions on Gamma.

The Neumann and Dirichlet boundary conditions can be functions of space (x and y) and time t.The periodic boundary conditions are used for the problem defined over a square of arbitrary side length L (the code is self contained and thus a grid does not need to be imported). If we label the sides of a square bn1, bn2, bn3 and bn4 then at each time step the boundary values on sides bn2 and bn1 are used to over-write the boundary values on sides bn4 and bn3 respectively, thus enforcing the periodicity of the boundary conditions.

There is no special format for entering the boundary functions, for example when running fe2dx_nd.m we might have:

       >> Enter the Dirichlet b.c. g1u(x,y,t) for u  0.2*sin(x-y)+t/10+0.2
       >> Enter the Dirichlet b.c. g1v(x,y,t) for v  0.3*cos(x-y)+t/20+0.3
       >> Enter the Neumann b.c. g2v(x,y,t) for v  0.0
       >> Enter the Neumann b.c. g2v(x,y,t) for v  0.6
    
## Some practical issues

Firstly, bear in mind that if you run a simulation with a large domain size and large final time T, coupled with small temporal and spatial discretization parameters, then the run-time in MATLAB can be prohibative.

Another point concerns the choice of parameters alpha, beta, and gamma used to run the code. In order for the local kinetics of the systems to be biologically meaningful, there are restrictions on the parameters that need to be satisfied (for further details see the references below).

Another issue concerns the stability and accuracy of the numerical solutions. The user needs to choose the time step sufficiently small to obtain converged solutions. This can easily be checked by reducing the time-step until the solutions corresponding to scheme 1 and scheme 2 (see the references below) are the same. Experience has shown that with many examples the time step needs to be of the order of 1/384.

One final point concerns the code with the Robin boundary conditions. If you run the code with k1 or k2 large and positive on part of the boundary then this part of the boundary acts as a source for the domain. Thus solvers may fail due to growing solutions, which leads to ill-conditioned coefficient matrices. No such problems occur if we choose k1 or k2 large and negative (although solutions may go negative if the time step is too large).

## Nonconvergence of GMRES

The GMRES algorithm does not always converge, or it may converge very slowly. One approach that we have incorporated into the codes to help overcome this problem is to use MATLAB's 'incomplete lu factorization' algorithm ILU (which replaces the previously used LUINC function). This provides a preconditioner for GMRES. See MATLAB's Help page under GMRES for usage.

If this still doesn't help you can try other preconditioners, or use a less sophisticated iterative method for solving the linear systems. For example, the Jacobi and Gauss-Seidel algorithms are widely used, although the run-times will likely be considerably longer than for GMRES. As these methods are not part of the standard suite of MATLAB functions they are provided below. The algorithms require that all diagonal elements aii of the coefficient matrix A are non-zero, and they will not always converge even then.

Files you may copy include:

    jacobi.m,    MATLAB code for iteratively solving the (square) system of linear equations Ax = b.
    gauss_seidel.m,    A variation (usually an improvement) of JACOBI.

## References

    Garvie M.R. , "Finite Difference Schemes for Reaction-Diffusion Equations Modeling Predator-Prey Interactions in MATLAB," Bulletin of Mathematical Biology (2007) Vol.69, No.3., pp. 931-956

    Garvie M.R. , Burkardt J., Morgan J. : "Simple Finite Element Methods for Approximating Predator-Prey Dynamics in Two Dimensions using MATLAB," Bulletin of Mathematical Biology (2015) Vol.77, No.3, pp. 548-578.

## Download the 'slow' codes for FE2D

These are the codes described in the 2nd paper referenced above. For more computationally intensive problems the 'fast' versions of these codes are recommended (see below). Files you may copy include:

    fe2d_n.m,  Scheme 2 applied to Kinetics 1 with pure Neumann boundary conditions.
    fe2dx_n.m,  Scheme 1 applied to Kinetics 1 with pure Neumann boundary conditions.
    fe2d_d.m,  Scheme 2 applied to kinetics 1 with pure Dirichlet boundary conditions.
    fe2dx_d.m,  Scheme 1 applied to kinetics 1 with pure Dirichlet boundary conditions.
    fe2d_r.m,  Scheme 2 applied to kinetics 1 with pure Robin boundary conditions.
    fe2dx_r.m,  Scheme 1 applied to kinetics 1 with pure Robin boundary conditions.
    fe2d_p.m,  Scheme 2 applied to kinetics 1 with periodic boundary conditions over the square.
    fe2dx_p.m,  Scheme 1 applied to kinetics 1 with periodic boundary conditions over the square.
    fe2d_nd.m,  Scheme 2 applied to kinetics 1 with mixed Neumann/Dirichlet boundary conditions.
    fe2dx_nd.m,  Scheme 1 applied to kinetics 1 with mixed Neumann/Dirichlet boundary conditions.
    fe2d_nr.m,  Scheme 2 applied to kinetics 1 with mixed Neumann/Robin boundary conditions.
    fe2dx_nr.m,  Scheme 1 applied to kinetics 1 with mixed Neumann/Robin boundary conditions.
    fe2dx_nr_alt.m,  Same as fe2dx_nr.m, but using an implicit approximation of the boundary terms
    mycopyright.txt,    FE2D copyright details.

## Download the 'fast' codes for FE2D

These are the same as the codes listed above, but optimized for speed. Files you may copy include the MATLAB M-files, and front ends with test data for running the fast codes. See the comments in the fast M-files for further implementation details:

    fe2d_n_fast.m,  Scheme 2 applied to Kinetics 1 with pure Neumann boundary conditions.
    fe2d_n_fast_test.m,  A front end with some test data that can be used to run fe2d_n_fast.m.
    fe2dx_n_fast.m,  Scheme 1 applied to Kinetics 1 with pure Neumann boundary conditions.
    fe2dx_n_fast_test.m,  A front end with some test data that can be used to run fe2dx_n_fast.m.
    fe2d_d_fast.m,  Scheme 2 applied to kinetics 1 with pure Dirichlet boundary conditions.
    fe2d_d_fast_test.m,  A front end with some test data that can be used to run fe2d_d_fast.m.
    fe2dx_d_fast.m,  Scheme 1 applied to kinetics 1 with pure Dirichlet boundary conditions.
    fe2dx_d_fast_test.m,  A front end with some test data that can be used to run fe2dx_d_fast.m.
    fe2d_r_fast.m,  Scheme 2 applied to kinetics 1 with pure Robin boundary conditions.
    fe2d_r_fast_test.m,  A front end with some test data that can be used to run fe2d_r_fast.m.
    fe2dx_r_fast.m,  Scheme 1 applied to kinetics 1 with pure Robin boundary conditions.
    fe2dx_r_fast_test.m,  A front end with some test data that can be used to run fe2dx_r_fast.m.
    fe2d_p_fast.m,  fe2d_p_fast.pdf    Scheme 2 applied to kinetics 1 with Periodic boundary conditions.
    fe2d_p_fast_test.m,  A front end with some test data that can be used to run fe2d_p_fast.m.
    fe2dx_p_fast.m,  Scheme 1 applied to kinetics 1 with Periodic boundary conditions.
    fe2dx_p_fast_test.m,  A front end with some test data that can be used to run fe2dx_p_fast.m.
    fe2d_nd_fast.m,  Scheme 2 applied to kinetics 1 with mixed Neumann/Dirichlet boundary conditions.
    fe2d_nd_fast_test.m,  A front end with some test data that can be used to run fe2d_nd_fast.m.
    fe2dx_nd_fast.m,  Scheme 1 applied to kinetics 1 with mixed Neumann/Dirichlet boundary conditions.
    fe2dx_nd_fast_test.m,  A front end with some test data that can be used to run fe2dx_nd_fast.m.
    fe2d_nr_fast.m,  Scheme 2 applied to kinetics 1 with mixed Neumann/Robin boundary conditions.
    fe2d_nr_fast_test.m,  A front end with some test data that can be used to run fe2d_nr_fast.m.
    fe2dx_nr_fast.m,  Scheme 1 applied to kinetics 1 with mixed Neumann/Robin boundary conditions.
    fe2dx_nr_fast_test.m,  A front end with some test data that can be used to run fe2dx_nr_fast.m.
    fe2dx_nr_alt_fast.m,  Same as fe2dx_nr_fast.m, with an implicit approximation of the boundary terms
    fe2dx_nr_alt_fast_test.m,  A front end with some test data that can be used to run fe2dx_nr_alt_fast.m.
    mycopyright.txt,    FE2D copyright details.

## Associated codes you will need

Files you may copy include:

    boundedges.m,    Finds all the boundary edges in a triangular mesh.
    subsetconnectivity.m    Finds the boundary edges in a triangular mesh for a portion of the boundary.
    timestamp.m,    Prints the current YMDHMS date as a timestamp.

PRED_PREY_ARB_FE2D is distributed under the GNU GPL; see the mycopyright.txt file for more information. 
