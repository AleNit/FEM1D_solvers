
% solve 1D Heat equation with the standard Galerkin method: 
%       du/dt + alp*d^2(u)/dz^2=f, for x in [a,b]
%       du(a)/dz=lbc;
%       u(b)=rbc;
% quadratic elements, Cranck-Nicholson time scheme

% A. Nitti (2024), Polytechnic University of Bari

clc
clear
close all


%% input parameters
a=0;                    % domain boundaries
b=1;
te=0.4;                 % final integration time
f = 0;                  % right hand side value
lbc = 0.0;             % boundary value (Neumann condition)
rbc = 0.0;              % boundary value (Dirichelet condition) 
ne=50;                   % number of elements
% nt=200;                 % numer of time steps
nt=600;                 % numer of time steps
u0=10.0;
alp=-1;


%% build stiffness matrix, mass matrix and right-hand-side
dom=[a,b];
dt=te/nt;
ndof=ne*2+1;                            % number of degrees-of- freedom
z = linspace(dom(1),dom(2),ndof)';      % nodal coordinates
[M,K,RHS]=getmats(z,2,ne,ndof,f);
K=alp.*K;


%% initialize arrays for time integration
un=u0.*ones(ndof,1);
un(ndof)=rbc;
t=0;


%% time integration
figure(1);
za=linspace(a,b,250);
narm=1000;  % number of harmonics for the analytical solution

for n=1:nt

    % compute analitical solution    
    ua=zeros(1,length(za));
    for j=1:narm
        ua=ua+(-1)^(j+1)/(2*j-1).*cos((2*j-1)/2*pi.*za).*exp(-(2*j-1)^2*pi^2*t/4);
    end
    ua=4/pi*u0.*ua;
    plot(za, ua, '-k'); hold on

    % plot current solution    
    plot(z, un, 'or');
    xlabel('z','interpreter','latex','fontsize',14);
    ylabel('u','interpreter','latex','fontsize',14);
    legend('analytical','numerical','interpreter','latex','fontsize',12, ...
        'location','southwest')
    str=strcat('t=',num2str(t));
    title(str,'interpreter','latex','fontsize',14)
    axis([dom(1),dom(2),rbc,u0+1])
    drawnow 
    hold off

    % assemble updated right hand side    
    RHSd=(M+0.5*dt.*K)*un+dt.*RHS;    
    Kd=M-0.5*dt.*K;

    % reduced stiffness matrix
    Kr=Kd(1:ndof-1,1:ndof-1);
    
    % compute reduced RHS
    Kbc=Kd(1:ndof-1,ndof);
    g=zeros(ndof-1,1);
    g(1)=lbc;
    RHSr=RHSd(1:ndof-1)-Kbc*rbc+dt.*alp.*g;
    
    % linear solver
    ur=linsolve(Kr,RHSr);
    
    % assemble full solution vector
    unp1=[ur; rbc];    

    un=unp1;
    t=t+dt;
    disp(['i=',num2str(n),', time=',num2str(t)])

end

