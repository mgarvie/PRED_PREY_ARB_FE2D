function fe2d_nd
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                               Discussion
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 'fe2d_nd.m'   2D finite element Matlab code for Scheme 2 applied
% to the predator-prey system with Kinetics 1. The nodes and elements
% of the unstructured grid are loaded from external files 't_triang.dat'
% and 'p_coord.dat' respectively, as are the list of nodes on which 
% Dirichlet and Neumann b.c.'s are to be imposed (from 'bn1_nodes.dat' and 
% 'bn2_nodes.dat' respectively).
%
% Boundary conditions:
%   Gamma1: Dirichlet
%   Gamma2: Neumann
%
% (C) 2009 Marcus R. Garvie. See 'mycopyright.txt' for details.
%
% Modified April 7, 2014
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                      Enter data for mesh geometry
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Read in 'p(2,n)', the 'n' coordinates of the nodes
load p_coord.dat -ascii
p = (p_coord)';
% Read in 't(3,no_elems)', the list of nodes for 'no_elems' elements
load t_triang.dat -ascii
t = (round(t_triang))';
% Read in 'bn1(1,isn1)', the nodes on Gamma1
load bn1_nodes.dat -ascii
bn1 = (round(bn1_nodes))';
% Read in 'bn2(1,isn2)', the nodes on Gamma2
load bn2_nodes.dat -ascii
bn2 = (round(bn2_nodes))';
% Construct the connectivity for the nodes on Gamma2
cpp = subsetconnectivity (t', p', bn2');
% Number of edges on Gamma2
[e2,junk] = size(cpp);
% Degrees of freedom per variable (n)
[junk,n]=size(p);
% Number of elements (no_elems)
[junk,no_elems]=size(t);
% Number of nodes on boundary Gamma1 (isn1)
[junk,isn1]=size(bn1);
% Extract vector of 'x' and 'y' values
x = p(1,:); y = p(2,:);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                     Enter data for model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% User inputs of parameters
alpha = input('Enter parameter alpha   ');
beta = input('Enter parameter beta   ');
gamma = input('Enter parameter gamma   ');
delta = input('Enter parameter delta   ');
T = input('Enter maximum time T   ');
delt = input('Enter time-step Delta t   ');
% User inputs of initial data
u0_str = input('Enter initial data function u0(x,y)   ','s');
u0_anon = @(x,y)eval(u0_str);   % create anonymous function
u = arrayfun(u0_anon,x,y)';
v0_str = input('Enter initial data function v0(x,y)   ','s');
v0_anon = @(x,y)eval(v0_str);   % create anonymous function
v = arrayfun(v0_anon,x,y)';
% Enter the boundary conditions
g1u_str = input('Enter the Dirichlet b.c. g1u(x,y,t) for u  ','s');
g1u = @(x,y,t)eval(g1u_str);   % create anonymous function
g1v_str = input('Enter the Dirichlet b.c. g1v(x,y,t) for v  ','s');
g1v = @(x,y,t)eval(g1v_str);   % create anonymous function
g2u_str = input('Enter the Neumann b.c. g2u(x,y,t) for u ','s');
g2u = @(x,y,t)eval(g2u_str);   % create anonymous function
g2v_str = input('Enter the Neumann b.c. g2v(x,y,t) for v ','s');
g2v = @(x,y,t)eval(g2v_str);   % create anonymous function
% Calculate and assign some constants                 
N=round(T/delt);
% Degrees of freedom per variable (n)
[junk,n]=size(p);
% Number of elements (no_elems)
[junk,no_elems]=size(t);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                               Assembly
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
m_hat=zeros(n,1);
K=sparse(n,n);
for elem = 1:no_elems
    % Identify nodes ni, nj and nk in element 'elem'
    ni = t(1,elem);
    nj = t(2,elem);
    nk = t(3,elem);
    % Identify coordinates of nodes ni, nj and nk
    xi = p(1,ni);
    xj = p(1,nj);
    xk = p(1,nk);
    yi = p(2,ni);
    yj = p(2,nj);
    yk = p(2,nk); 
    % Calculate the area of element 'elem'
    triangle_area = abs(xj*yk-xk*yj-xi*yk+xk*yi+xi*yj-xj*yi)/2;
    % Calculate some quantities needed to construct elements in K
    h1 = (xi-xj)*(yk-yj)-(xk-xj)*(yi-yj);
    h2 = (xj-xk)*(yi-yk)-(xi-xk)*(yj-yk);
    h3 = (xk-xi)*(yj-yi)-(xj-xi)*(yk-yi);
    s1 = (yj-yi)*(yk-yj)+(xi-xj)*(xj-xk);
    s2 = (yj-yi)*(yi-yk)+(xi-xj)*(xk-xi);
    s3 = (yk-yj)*(yi-yk)+(xj-xk)*(xk-xi);
    t1 = (yj-yi)^2+(xi-xj)^2;  % g* changed to t*
    t2 = (yk-yj)^2+(xj-xk)^2;
    t3 = (yi-yk)^2+(xk-xi)^2;
    % Calculate local contributions to m_hat 
    m_hat_i = triangle_area/3;
    m_hat_j = m_hat_i;
    m_hat_k = m_hat_i;
    % calculate local contributions to K
    K_ki = triangle_area*s1/(h3*h1);
    K_ik = K_ki;
    K_kj = triangle_area*s2/(h3*h2);
    K_jk = K_kj;
    K_kk = triangle_area*t1/(h3^2);
    K_ij = triangle_area*s3/(h1*h2);
    K_ji = K_ij;
    K_ii = triangle_area*t2/(h1^2);
    K_jj = triangle_area*t3/(h2^2);
    % Add contributions to vector m_hat
    m_hat(nk)=m_hat(nk)+m_hat_k;
    m_hat(nj)=m_hat(nj)+m_hat_j;
    m_hat(ni)=m_hat(ni)+m_hat_i;
    % Add contributions to K
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
% Construct matrix L
ivec=1:n;
IM_hat=sparse(ivec,ivec,1./m_hat,n,n);
L=delt*IM_hat*K;
% Construct matrices B1 & B2
B1=sparse(1:n,1:n,1,n,n)+L;
B2=sparse(1:n,1:n,1,n,n)+delta*L;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           Time-stepping procedure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for nt=1:N
    tn = nt*delt;
    % Evaluate modified functional response
    hhat = u./(alpha + abs(u));
    % Update right-hand-side of linear system 
    F = u - u.*abs(u) - v.*hhat;
    G = beta*v.*hhat - gamma*v;
    rhs_u = u + delt*F;
    rhs_v = v + delt*G;    
    % Impose Neumann boundary condition on Gamma2
    for i = 1:e2
        node1 = cpp(i,1);
        node2 = cpp(i,2);
        x1 = p(1,node1);
        y1 = p(2,node1);
        x2 = p(1,node2);
        y2 = p(2,node2);
        im_hat1 = 1/m_hat(node1);
        im_hat2 = 1/m_hat(node2);
        gamma12 = sqrt((x1-x2)^2 + (y1-y2)^2);
        rhs_u(node1) = rhs_u(node1) + delt*g2u(x1,y1,tn)*im_hat1*gamma12/2;
        rhs_u(node2) = rhs_u(node2) + delt*g2u(x2,y2,tn)*im_hat2*gamma12/2;
        rhs_v(node1) = rhs_v(node1) + delt*g2v(x1,y1,tn)*im_hat1*gamma12/2;
        rhs_v(node2) = rhs_v(node2) + delt*g2v(x2,y2,tn)*im_hat2*gamma12/2;
    end
    % Impose dirichlet boundary conditions on Gamma1
    for i = 1:isn1
        node = bn1(i);
        xx = p(1,node);
        yy = p(2,node);
        B1(node,:)=0; 
        B1(node,node)=1;
        rhs_u(node)=g1u(xx,yy,tn);
        B2(node,:)=0; 
        B2(node,node)=1;
        rhs_v(node)=g1v(xx,yy,tn);
    end
    % Do the LU factorization of B1 and B2
    [LB1,UB1] = ilu(B1,struct('type','ilutp','droptol',1e-5));
    [LB2,UB2] = ilu(B2,struct('type','ilutp','droptol',1e-5));
    % Solve for u and v using GMRES
    [u,flagu,relresu,iteru]=gmres(B1,rhs_u,[],1e-6,[],LB1,UB1,u);
    if flagu~=0 flagu,relresu,iteru,error('GMRES did not converge'),end
    [v,flagv,relresv,iterv]=gmres(B2,rhs_v,[],1e-6,[],LB2,UB2,v);
    if flagv~=0 flagv,relresv,iterv,error('GMRES did not converge'),end   
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                               Plot solutions 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot solution for u
figure;
set(gcf,'Renderer','zbuffer');
trisurf(t',x,y,u,'FaceColor','interp','EdgeColor','interp');
colorbar;axis off;title('u');
view ( 2 );
axis equal on tight;
% Plot solution for v
figure;
set(gcf,'Renderer','zbuffer');
trisurf(t',x,y,v,'FaceColor','interp','EdgeColor','interp');
colorbar;axis off;title('v');
view ( 2 );
axis equal on tight;