D_LSCF = transp([9.58E-18;9.96E-17;7.64E-16;4.58E-15;2.24E-14;9.20E-14;3.27E-13;
    1.03E-12; 2.90E-12]);
D_GDC = transp([1.59E-14; 4.59E-14; 1.14E-13; 2.61E-13; 5E-13; 1E-12; 2E-12; 3E-12;
    5E-12]);
k_LSCF = transp([8.40E-11; 1.76E-10; 3.54E-10;  6.05E-10; 1.17E-09; 4.69E-09; 1.42E-08
    4.45E-08; 1.14E-07]);
for n = 1:9
D1=D_LSCF(n); % Define diffusivity of cathode material [um^2/s] (LSCF)
D2=D_GDC(n);  % Define diffusivity of electrolyte material [um^2/s] (GDC)
k = k_LSCF(n); 
r = 0; %Define interfacial resistance [s/m]
t_in = 10.5;
ProfileLength = 200;
InterfacePosition = 25;

[X,pro]=InterfaceCNop_inline(D1,D2,k,r,t_in,ProfileLength,InterfacePosition);

plot (X*1e6,pro)
hold on
end

function [X,pro]=InterfaceCNop_inline(D1,D2,k,r,t_in,ProfileLength,InterfacePosition)
PixelNo=139;
dx=ProfileLength*1e-6/(PixelNo-1); %spatial step
X=0:dx:ProfileLength*1e-6; %domain
t=t_in*3600; %s %initial exchange duration
h=k/D1;
C_bg = 0;
C_gas = 1;
% Crank1=(C_bg+(C_gas-C_bg)*...
%     (erfc(X./(2*sqrt(D1*t)))-...
%     exp(h.*X+t*D1*h^2)...
%     .*erfc(X./(2*sqrt(D1*t))+h*sqrt(D1*t))))';

% dt=5*round(dx^2/(2*min(D1,D2)),3,'significant'); % Time step
dt = 1;
Nt=t/dt+1; % Number of time steps
Nx=ProfileLength*1e-6/dx+1; % Define number of spatial nodes

L = 1e-6*[InterfacePosition;ProfileLength-InterfacePosition]; % Vector of the layer lengths, first element is the first layer
D_int=2/(2*r/dx+1/D1+1/D2);

sigma_1=D1*(dt)/(2*dx^2); % Calculate sigma of LSCF
sigma_2=D2*(dt)/(2*dx^2); % Calculate sigma of GDC
sigma_int=D_int*(dt)/(2*dx^2);% Calculate sigma of the first interface


%%

% Time stepping
%% Build simulation components
% Set up a sub-matrix L which contains the information regarding the surface exchange as well as diffusion through the LSCF layer

% Set up the vector for the subdiagonal
sub = zeros(round(Nx)-1,1);
sub(1:round(L(1)/dx)-1,1)=sigma_1;
sub(round(L(1)/dx):round(L(1)/dx)+1,1)=sigma_int;
sub(round(L(1)/dx)+1:end,1)=sigma_2;

% Set up the vector for the superdiagonal
sup = zeros(round(Nx)-1,1);
sup(1:round(L(1)/dx)-1,1)=sigma_1;
sup(round(L(1)/dx)-0:round(L(1)/dx),1)=sigma_int;
sup(round(L(1)/dx)+1:end,1)=sigma_2;
% Set up the vector for the central diagonal
s1 = zeros(round(Nx),1);
s1(1:round(L(1)/dx)-1,1)=1-2*sigma_1;
s1(round(L(1)/dx),1)=1-(sigma_int+sigma_1);
s1(round(L(1)/dx)+1,1)=1-(sigma_int+sigma_2);
s1(round(L(1)/dx)+2:end,1)=1-2*sigma_2;

s1_n = zeros(round(Nx),1);
s1_n(1:round(L(1)/dx)-1,1)=1+2*sigma_1;
s1_n(round(L(1)/dx),1)=1+(sigma_int+sigma_1);
s1_n(round(L(1)/dx)+1,1)=1+(sigma_int+sigma_2);
s1_n(round(L(1)/dx)+2:end,1)=1+2*sigma_2;

% Create the matrix at the current timestep (Mirror boundary condition)
A = full(gallery('tridiag',sub,s1,sup));
A(1,1:2) = [1-2*sigma_1*(1+dx*h) 2*sigma_1]; %Surface boundary condition
A(end,end-2:end) = [sigma_2 -2*sigma_2 1+sigma_2]; %Mirror boundary condition parameters

% Create the matrix at the future timestep
A_n = full(gallery('tridiag',-sub,s1_n,-sup));
A_n(1,1:2) = [1+2*sigma_1*(1+dx*h) -2*sigma_1]; %Surface boundary condition
A_n(end,end-2:end) = [-sigma_2 2*sigma_2 1-sigma_2]; %Mirror boundary condition parameters

% Create the matrix at the current timestep (Dirichlet boundary condition)
A_Di = full(gallery('tridiag',sub,s1,sup));
A_Di(1,1:2) = [1-2*sigma_1*(1+dx*h), 2*sigma_1]; %Surface boundary condition
A_Di(end,end-2:end) = [0 0 1]; % Dirichlet boundary condition

A_n_Di = full(gallery('tridiag',-sub,-s1+2,-sup));
A_n_Di(1,1:2) = [1+2*sigma_1*(1+dx*h), -2*sigma_1]; %Surface boundary condition
A_n_Di(end,end-2:end) = [0 0 1]; %Mirror boundary condition parameters

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Define surface exchange vector
G = zeros(round(Nx), 1);
G(1) = 4*dx*sigma_1*h; %surface boundary condition vector G
C = zeros(round(Nx),1);%Initial position vector, 0 everywhere
C_Di = zeros(round(Nx),1);% Initial position vector for the Dirichlet bc


%% Run the iterative matrix calculation

A_nI = inv(A_n);
A_nI_Di = inv(A_n_Di);

for t=1:Nt
    C = A_nI*(A*C+G);
end

if C(end) < 0.001
    pro = C;
else
    
    for t=1:Nt
        C_Di = A_nI_Di*(A_Di*C_Di+G);
    end
    
    pro = (C+C_Di)/2;
end


end