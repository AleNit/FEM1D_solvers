
% solve 1D Poisson equation with the standard Galerkin method: 
%       d^2(u)/dx^2=f, for x in [a,b]
%       du(a)/dx=lbc;
%       u(b)=rbc;

% A. Nitti (2024), Polytechnic University of Bari

clc
clear
close all


%% input parameters
a=0;                    % domain boundaries
b=2.5;
f = 1;                  % right hand side value
lbc = -0.5;             % boundary value (Neumann condition)
rbc = 1.0;              % boundary value (Dirichelet condition) 
ne=5;                   % number of elements


%% build stiffness matrix, mass matrix and right-hand-side
dom=[a,b];
ndof=ne*2+1;                            % number of degrees-of-freedom
x = linspace(dom(1),dom(2),ndof)';      % nodal coordinates
[~,K,RHS]=getmats(x,2,ne,ndof,f);
K=-K;                                   % due to integration by parts


%% apply boundary conditions
Kr=K(1:ndof-1,1:ndof-1);                      % reduced stiffness matrix

Kbc=K(1:ndof-1,ndof);
g=zeros(ndof-1,1);
g(1)=lbc;
RHSr=RHS(1:ndof-1)-Kbc*rbc+g;              % reduced RHS

% linear solver
ur=linsolve(Kr,RHSr);

% assemble full solution vector
u=[ur; rbc];


%% analytical solution
xa=linspace(dom(1),dom(2),200);
ua=f/2.*(dom(2)-xa).*(2*dom(1)-dom(2)-xa)+lbc.*(xa-dom(2))+rbc;


%% plot solution
figure(1);
plot(x, u, 'or');
hold on
plot(xa, ua, '-k');
xlabel('x','interpreter','latex','fontsize',14);
ylabel('u','interpreter','latex','fontsize',14);
legend('numerical','analytical','interpreter','latex','fontsize',12, ...
    'location','southeast')
xlim([dom(1),dom(2)])


%% compute L^2-norm error vs the manufactured analytical solution
ua2=f/2.*(dom(2)-x).*(2*dom(1)-dom(2)-x)+lbc.*(x-dom(2))+rbc;
error = sqrt(sum((ua2-u).^2)/sum(ua2.^2));

disp(['L2 error norm: ',num2str(error)])

