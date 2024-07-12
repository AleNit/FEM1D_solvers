
% build stiffness matrix, mass matrix and right-hand-side
% quadratic elements with Gauss-Legendre quadrature are used

function [M,K,RHS]=getmats(coord,p,ne,n,F)

% create node connectivity
In=[(1:p:n-p)',(2:p:n-p+1)',(3:p:n)'];

% Gauss point location and weights in the parent space [-1,1]
gp=[-sqrt(3/5),0,sqrt(3/5)];
gw=[5/9,8/9,5/9];

% define shape function values and derivatives at quadrature nodes 
% (2nd degree Lagrange polynomials in the parent space)
sf=zeros(p+1,p+1);
sfd=sf;
for i = 1:p+1
    sf(1,i)  = 0.5*gp(i)^2+0.5*gp(i);
    sf(2,i)  = -gp(i)^2+1;
    sf(3,i)  = 0.5*gp(i)^2-0.5*gp(i);
    sfd(1,i) = gp(i)+0.5;
    sfd(2,i) = -2*gp(i);
    sfd(3,i) = gp(i)-0.5;
end

% initialize global matrices and right-hand-side vector
K = zeros(n);
M=K;
RHS = zeros(n, 1);

% element loop
for i=1:ne

    Ine=In(i,:);
    zn=coord(Ine);
    h=zn(end)-zn(1);
    J=h/2;
    
    % element matrices
    Ke=zeros(p+1,p+1);
    Me=Ke;
    RHSe=zeros(p+1,1);
    for j=1:p+1 
        phi=sf(:,j);
        dphidxi=sfd(:,j);    
        Ke   = Ke+dphidxi*dphidxi'./J*gw(j);
        Me   = Me+phi*phi'.*J*gw(j);
        RHSe = RHSe+phi*F.*J*gw(j);
    end

    % plug into global stiffness matrix
    K(Ine,Ine) = K(Ine,Ine)+Ke;
    M(Ine,Ine) = M(Ine,Ine)+Me;
    RHS(Ine)   = RHS(Ine)+RHSe;
    
end

end
