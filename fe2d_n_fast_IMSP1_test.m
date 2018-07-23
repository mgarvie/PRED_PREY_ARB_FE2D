function fe2d_n_fast_IMSP1_test ( )

%*****************************************************************************80
%
%% FE2D_N_FAST_TEST_IMSP1 tests the FE2D_N_FAST_IMSP1 code.
%
%  Discussion:
%
%    This function sets all parameter values and initial condition information
%    necessary to execute the "fast" version of the fe2d_n algorithm.
%
%  Licensing:
%
%    Copyright (C) 2014 Marcus R. Garvie. 
%    See 'mycopyright.txt' for details.
%
%  Modified:
%
%    28 April 2014
%
%  Authors:
%
%    Marcus R. Garvie and John Burkardt. 
%
%  Reference:
%
%    Marcus R Garvie, John Burkardt, Jeff Morgan,
%    Simple Finite Element Methods for Approximating Predator-Prey Dynamics
%    in Two Dimensions using MATLAB,
%    Bulletin of Mathematical Biology (2015) Vol.77, No.3, pp. 548-578. 
%
%  Modified 19 June 2016
%
%  Authors:
%
%   Fasma Diele and Marcus R. Garvie.
%
%  Reference:
%
%    Fasma Diele, Marcus R. Garvie, Catalin Trenchea
%    Analysis of first order in time implicit-symplectic 
%    scheme for predator-prey systems
%
%    Submitted, 2016.
%
%
  timestamp ( );
  fprintf ( 1, '\n' );
  fprintf ( 1, 'FE2D_N_FAST_IMSP1_TEST:\n' );
  fprintf ( 1, '  Test the FE2D_N_FAST_IMSP1 function\n' );
  fprintf ( 1, '  which applies Neumann boundary conditions as it\n' );
  fprintf ( 1, '  approximates a solution to a predator-prey system.\n' );
%
%  Set the parameters.
%
  alpha = 0.4;
  beta = 2.0;
  gamma = 0.6;
  Du = 1.0;
  Dv = 1.0;
  
  T = 150;
  delt = 1.0 / 3.0;
%    delt = 1.0 / 24.0;
%    delt = 1.0 / 384.0;

  t = tic;
  [u,v]=fe2d_n_fast_IMSP1 ( alpha, beta, gamma, Du, Dv, T, delt, @u0f, @v0f, @guf, @gvf );
  t = toc ( t );



  fprintf ( 1, '  Execution took %10.2g minutes \n', t / 60.0 );
%
%  Terminate.

  fprintf ( 1, '\n' );
  fprintf ( 1, 'FE2D_N_FAST_TEST_IMSP1:\n' );
  fprintf ( 1, '  Normal end of execution.\n' );
  fprintf ( 1, '\n' );
  timestamp ( );

  return
end

function value = u0f ( x, y )

%*****************************************************************************80
%
%% U0F evaluates the initial condition for U.
%
%  Licensing:
%
%    Copyright (C) 2014 Marcus R. Garvie. 
%    See 'mycopyright.txt' for details.
%
%  Modified:
%
%    26 April 2014
%
%  Author:
%
%    Marcus R. Garvie. 
%
%  Parameters:
%
%    Input, real X, Y, a location in the region.
%
%    Output, real VALUE, the initial condition for U at (X,Y).
%
  value = 6.0 / 35.0 - 2.0E-07 * ( x - 0.1 * y - 225.0 ) * ( x - 0.1 * y - 675.0 );

  return
end

function value = v0f ( x, y )

%*****************************************************************************80
%
%% V0F evaluates the initial condition for V.
%
%  Licensing:
%
%    Copyright (C) 2014 Marcus R. Garvie. 
%    See 'mycopyright.txt' for details.
%
%  Modified:
%
%    26 April 2014
%
%  Author:
%
%    Marcus R. Garvie. 
%
%  Parameters:
%
%    Input, real X, Y, a location in the region.
%
%    Output, real VALUE, the initial condition for V at (X,Y).
%
  value = 116.0 / 245.0 - 3.0E-05 * ( x - 450.0 ) - 1.2E-04 * ( y - 150.0 );

  return
end

function value = guf ( x, y, t )

%*****************************************************************************80
%
%% GUF evaluates the Neumann boundary condition for U.
%
%  Licensing:
%
%    Copyright (C) 2014 Marcus R. Garvie. 
%    See 'mycopyright.txt' for details.
%
%  Modified:
%
%    28 April 2014
%
%  Author:
%
%    Marcus R. Garvie. 
%
%  Parameters:
%
%    Input, real X, Y, a location on the boundary.
%
%    Input, real T, the time.
%
%    Output, real VALUE, the prescribed value of dU/dn at (X,Y,T).
%
  value = 0.0;

  return
end
function value = gvf ( x, y, t )

%*****************************************************************************80
%
%% GVF evaluates the Neumann boundary condition for V.
%
%  Licensing:
%
%    Copyright (C) 2014 Marcus R. Garvie. 
%    See 'mycopyright.txt' for details.
%
%  Modified:
%
%    28 April 2014
%
%  Author:
%
%    Marcus R. Garvie. 
%
%  Parameters:
%
%    Input, real X, Y, a location on the boundary.
%
%    Input, real T, the time.
%
%    Output, real VALUE, the prescribed value of dV/dn at (X,Y,T).
%
  value = 0.0;

  return
end
