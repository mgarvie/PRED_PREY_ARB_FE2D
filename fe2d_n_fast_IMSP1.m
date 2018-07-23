function [u,v]=fe2d_n_fast_IMSP1 ( alpha, beta, gamma, Du, Dv, T, delt, u0f, v0f, guf, gvf )

%*****************************************************************************80
%
%% FE2D_N_FAST_IMSP1 applies IMSP first order with Kinetics 1 to predator prey in a region.
%
%  Discussion:
%
%    Neumann boundary conditions are applied.
%
%    The parameters ALPHA, BETA, GAMMA, Du and Dv appear in the predator-prey
%    equations as follows:
%
%      dUdT = Du * nabla U +      U*V/(U+ALPHA) + U*(1-U)
%      dVdT = Dv * nabla V + BETA*U*V/(U+ALPHA) - GAMMA * V
%
%  Licensing:
%
%    Copyright (C) 2014 Marcus R. Garvie. 
%
%  Modified:
%
%    29 April 2014
%
%  Authors:
%
%    Marcus R. Garvie and John Burkardt.
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

%*****************************************************************************80
%  Enter data for mesh geometry.
%*****************************************************************************80
%
%  Read in 'p(2,n)', the 'n' coordinates of the nodes.
%
  load p_coord.dat -ascii
  p = ( p_coord )';
%
%  Read in 't(3,no_elems)', the list of nodes for 'no_elems' elements.
%
  load t_triang.dat -ascii
  t = ( round ( t_triang ) )';
%
%  Construct the connectivity for the nodes on Gamma.
%
  edges = boundedges ( p', t' );
%
%  E = number of edges on Gamma.
%
  [ e, ~ ] = size ( edges );
%
%  N = degrees of freedom per variable.
%
  [ ~, n ] = size ( p );
%
%  NO_ELEMS = number of elements.
%
  [ ~, no_elems ] = size ( t );
%
%  Extract vector of 'x' and 'y' values.
%
  x = p(1,:); 
  y = p(2,:);

%*****************************************************************************80
%  Enter data for model.
%*****************************************************************************80

  if ( nargin < 1 )
    alpha = input ( 'Enter parameter alpha:  ' );
  elseif ( ischar ( alpha ) )
    alpha = str2num ( alpha );
  end

  if ( nargin < 2 )
    beta = input ( 'Enter parameter beta:  ' );
  elseif ( ischar ( beta ) )
    beta = str2num ( beta );
  end

  if ( nargin < 3 )
    gamma = input ( 'Enter parameter gamma:  ' );
  elseif ( ischar ( gamma ) )
    gamma = str2num ( gamma );
  end

  if ( nargin < 4 )
    Du = input ( 'Enter parameter Du:  ' );
  elseif ( ischar ( Du ) )
    Du = str2num ( Du );
  end
  
  if ( nargin < 5 )
    Dv = input ( 'Enter parameter Dv:  ' );
  elseif ( ischar ( Dv ) )
    Dv = str2num ( Dv );
  end
   
  if ( nargin < 6 )
    T = input ( 'Enter maximum time T:  ' );
  elseif ( ischar ( T ) )
    T = str2num ( T );
  end

  if ( nargin < 7 )
    delt = input ( 'Enter time-step delt:  ' );
  elseif ( ischar ( delt ) )
    delt = str2num ( delt );
  end
  
  

  fprintf ( 1, '  Using ALPHA = %g\n', alpha );
  fprintf ( 1, '  Using BETA = %g\n', beta );
  fprintf ( 1, '  Using GAMMA = %g\n', gamma );
  fprintf ( 1, '  Using Du = %g\n', Du );
  fprintf ( 1, '  Using Dv = %g\n', Dv );
  fprintf ( 1, '  Using T = %g\n', T );
  fprintf ( 1, '  Using DELT = %g\n', delt );
%
%  Initial conditions.
%
  if ( nargin < 8 )
    u0_str = input ( 'Enter initial data function u0(x,y):  ', 's' );
    u0f = @(x,y) eval ( u0_str );
  elseif ( ischar ( u0f ) )
    u0_str = u0f;
    u0f = @(x,y) eval ( u0_str );
  end

  u = ( arrayfun ( u0f, x, y ) )';

  if ( nargin < 9 )
    v0_str = input ( 'Enter initial data function v0(x,y):  ', 's' );
    v0f = @(x,y) eval ( v0_str );
  elseif ( ischar ( v0f ) )
    v0_str = v0f;
    v0f = @(x,y) eval ( v0_str );
  end

  v = ( arrayfun ( v0f, x, y ) )';
%
%  Boundary conditions.
%
  if ( nargin < 10 )
    gu_str = input('Enter the Neumann b.c. gu(x,y,t) for u  ','s');
    guf = @(x,y,t) eval(gu_str); 
  elseif ( ischar ( guf ) )
    gu_str = guf;
    guf = @(x,y,t) eval(gu_str); 
  end

  if ( nargin < 11 )
    gv_str = input('Enter the Neumann b.c. gv(x,y,t) for v  ','s');
    gvf = @(x,y,t)eval(gv_str);
  elseif ( ischar ( gvf ) )
    gv_str = gvf;
    gvf = @(x,y,t)eval(gv_str);
  end
%
%  N = number of time steps.
%
  N = round ( T / delt );
  fprintf ( 1, '  Taking N = %d time steps\n', N );

%*****************************************************************************80
%  Assembly.
%*****************************************************************************80

  m_hat = zeros ( n, 1 );
  K = sparse ( n, n );

  for elem = 1 : no_elems
%
%  Identify nodes ni, nj and nk in element 'elem'.
%
    ni = t(1,elem);
    nj = t(2,elem);
    nk = t(3,elem);
%
%  Identify coordinates of nodes ni, nj and nk.
%
    xi = p(1,ni);
    xj = p(1,nj);
    xk = p(1,nk);
    yi = p(2,ni);
    yj = p(2,nj);
    yk = p(2,nk);
%
%  Calculate the area of element 'elem'.
%
    triangle_area = abs(xj*yk-xk*yj-xi*yk+xk*yi+xi*yj-xj*yi)/2;
%
%  Calculate some quantities needed to construct elements in K.
%
    h1 = (xi-xj)*(yk-yj)-(xk-xj)*(yi-yj);
    h2 = (xj-xk)*(yi-yk)-(xi-xk)*(yj-yk);
    h3 = (xk-xi)*(yj-yi)-(xj-xi)*(yk-yi);
    s1 = (yj-yi)*(yk-yj)+(xi-xj)*(xj-xk);
    s2 = (yj-yi)*(yi-yk)+(xi-xj)*(xk-xi);
    s3 = (yk-yj)*(yi-yk)+(xj-xk)*(xk-xi);
    t1 = (yj-yi)^2+(xi-xj)^2;
    t2 = (yk-yj)^2+(xj-xk)^2;
    t3 = (yi-yk)^2+(xk-xi)^2;
%
%  Calculate local contributions to m_hat.
%
    m_hat_i = triangle_area/3;
    m_hat_j = m_hat_i;
    m_hat_k = m_hat_i;
%
%  Calculate local contributions to K.
%
    K_ki = triangle_area*s1/(h3*h1);
    K_ik = K_ki;
    K_kj = triangle_area*s2/(h3*h2);
    K_jk = K_kj;
    K_kk = triangle_area*t1/(h3^2);
    K_ij = triangle_area*s3/(h1*h2);
    K_ji = K_ij;
    K_ii = triangle_area*t2/(h1^2);
    K_jj = triangle_area*t3/(h2^2);
%
%  Add contributions to vector m_hat.
%
    m_hat(nk)=m_hat(nk)+m_hat_k;
    m_hat(nj)=m_hat(nj)+m_hat_j;
    m_hat(ni)=m_hat(ni)+m_hat_i;
%
%  Add contributions to K.
%
    K=K+sparse(nk,ni,K_ki,n,n);
    K=K+sparse(ni,nk,K_ik,n,n);
    K=K+sparse(nk,nj,K_kj,n,n);
    K=K+sparse(nj,nk,K_jk,n,n);
    K=K+sparse(nk,nk,K_kk,n,n);
    K=K+sparse(ni,nj,K_ij,n,n);
    K=K+sparse(nj,ni,K_ji,n,n);
    K=K+sparse(ni,ni,K_ii,n,n);
    K=K+sparse(nj,nj,K_jj,n,n); 
  end
%
%  Construct matrix L.
%
  ivec = 1 : n;
  IM_hat = sparse(ivec,ivec,1./m_hat,n,n);
  L = delt * IM_hat * K;
%
%  Construct matrices B1 and B2.
%
  B1 = sparse(1:n,1:n,1,n,n) + Du* L;
  B2 = sparse(1:n,1:n,1,n,n) + Dv * L;
%
%  Do the incomplete LU factorisation of B1 and B2.
%
  [ LB1, UB1 ] = ilu ( B1, struct('type','ilutp','droptol',1e-5) );
  [ LB2, UB2 ] = ilu ( B2, struct('type','ilutp','droptol',1e-5) );

%*****************************************************************************80
%  Time-stepping.
%*****************************************************************************80

  for nt = 1 : N

    tn = nt * delt;
%
%  Evaluate functional response.
%
    hhat = u ./ ( alpha +  u  );
%
%  Update right-hand-side of linear system.
%
    G = beta * hhat - gamma ;
    rhs_v = v./(1 - delt * G);
    F = u - u .*  u  - rhs_v .* hhat;
    rhs_u = u + delt * F;

%
%  Impose Neumann boundary condition on Gamma.
%
    for i = 1 : e
      node1 = edges(i,1);
      node2 = edges(i,2);
      x1 = p(1,node1);
      y1 = p(2,node1);
      x2 = p(1,node2);
      y2 = p(2,node2);
      im_hat1 = 1/m_hat(node1);
      im_hat2 = 1/m_hat(node2);
      gamma12 = sqrt((x1-x2)^2 + (y1-y2)^2);
      rhs_u(node1) = rhs_u(node1) + Du * guf (x1,y1,tn) * im_hat1*gamma12/2;
      rhs_u(node2) = rhs_u(node2) + Du * guf (x2,y2,tn) * im_hat2*gamma12/2;
      rhs_v(node1) = rhs_v(node1) + Dv * gvf (x1,y1,tn) * im_hat1*gamma12/2;
      rhs_v(node2) = rhs_v(node2) + Dv * gvf (x2,y2,tn) * im_hat2*gamma12/2;
    end
%
%  Solve for u and v using GMRES.
%
    [u,flagu,relresu,iteru] = gmres ( B1,rhs_u,[],1e-6,[],LB1,UB1,u );

    if flagu ~= 0 
      flagu
      relresu
      iteru
      error('GMRES did not converge')
    end

    [v,flagv,relresv,iterv] = gmres ( B2,rhs_v,[],1e-6,[],LB2,UB2,v );

    if flagv ~= 0 
      flagv
      relresv
      iterv
      error('GMRES did not converge')
    end   

  end

%*****************************************************************************80
%  Plot solutions.
%*****************************************************************************80
%
%  Plot U;
%
  figure;
  set(gcf,'Renderer','zbuffer');
  trisurf(t',x,y,u,'FaceColor','interp','EdgeColor','interp');
  colorbar;
  axis off;
  title('u');
  view ( 2 );
  axis equal on tight; 
  filename = 'fe2d_n_fast_u.png';
  print ( '-dpng', filename );
  fprintf ( 1, '  Saved graphics file "%s"\n', filename );
%
%  Plot V.
%
  figure;
  set(gcf,'Renderer','zbuffer');
  trisurf(t',x,y,v,'FaceColor','interp','EdgeColor','interp');
  colorbar;
  axis off;
  title('v');
  view ( 2 );
  axis equal on tight;
  filename = 'fe2d_n_fast_v.png';
  fprintf ( 1, '  Saved graphics file "%s"\n', filename );
  print ( '-dpng', filename );

  return
end
