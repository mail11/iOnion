clear
InterfacePosition = 25;
ProfileLength = 1000;
t_in = 0.5:0.5:24;
T = 400:100:800;
kB = 8.6173e-05;
Tinv = 1./(T+273.5);
v =  linspace(1,length(T),length(T))./linspace(400,800,length(T));
D_LSCF = 0.0049255*exp((-1.9659/kB)*Tinv);
D_GDC = 8.0252e-08*exp((-0.8968/kB)*Tinv);
k_LSCF = transpose(zeros(round(length(T)), 1));
k_LSCF(1:(end+1)/2) = 7.0107e-06*exp((-0.632/kB)*Tinv(1:(end+1)/2));
k_LSCF((end+1)/2+1:end) = 53.9649*exp((-1.772/kB)*Tinv((end+1)/2+1:end));
r = 0;
for t_idx = 1:length(t_in)
    for T_idx = 1:length(T)
        D1 = D_LSCF(round(v(T_idx)*T(T_idx)));
        D2 = D_LSCF(round(v(T_idx)*T(T_idx)));
        k = k_LSCF(round(v(T_idx)*T(T_idx)));
        [X,pro] = InterfaceCNop_inline(D1,D2,k,r,t_in(t_idx),ProfileLength,InterfacePosition);
        CPrimeInt(T_idx,t_idx) = pro(1);
    end
end

contour(t_in,T,CPrimeInt)
c = colorbar;
c.Label.String = 'C at the surface at r = inf with D1 = D2';
set(gca,'YDir','normal')
xlabel('t(h)')
ylabel('T(C)')


function [X,pro]=InterfaceCNop_inline(D1,D2,k,r,t_in,ProfileLength,InterfacePosition)
PixelNo=200;
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
% 
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