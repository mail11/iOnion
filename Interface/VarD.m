function C = VarD(Dk)

global data;
global L;
global delta_x;
global delta_t;
global steps_x;
global Duration;
global Layers;
global C_gas;
global CM;
global xq;

%% Define parameters
% Material properties
%D_1=1.0E-02; % Define diffusivity of cathode material [um^2/s] (LSCF)
D_2=1.0E-03;
%D_3=0.00001E-03;% Define diffusivity of electrolyte material [um^2/s] (GDC)
%k=5.8E-02; % Define surface exchange coefficient of LSCF [um/s]
r1 = 0;
r2 = 0; %Define interfacial resistance [s/um]
int_width=0;% Define the width of the interface region (0 is the default value)


% Experimental setup
Duration=1.667; % Time of exchange in hours (found with Kiloran correction)
L = [0.1;0.0;0.107]; % Vector of the layer lengths, first element is the first layer

% Simulation parameters
delta_x=0.001; % Define the spatial step [um]
delta_t=2; % Define the time step [s] 
x_s = 0:delta_x*10e2:10;
f = 1-((1+exp(-Dk(4)))^-1)*exp(-x_s/(exp(Dk(5))));
Dv = Dk(1)*f;
sigma_v = Dv*(delta_t)/(2*delta_x^2);

% Derived properties
D_int1=2/(2*r1/delta_x+1/Dk(1)+1/D_2); 
D_int2=2/(2*r2/delta_x+1/D_2+1/Dk(2));
Length=sum(L); % Total sample length [um]
steps_x=Length/delta_x+1; % Define number of spatial nodes
steps_t=Duration*3600/delta_t+1; % Define number of time steps
x=0:delta_x:Length; % Define the length of the profile in steps of delta_x
sigma_1=Dk(1)*(delta_t)/(2*delta_x^2); % Calculate sigma of LSCF
sigma_2=D_2*(delta_t)/(2*delta_x^2); % Calculate sigma of GDC
sigma_3=Dk(2)*(delta_t)/(2*delta_x^2); % Calculate sigma of GDC
sigma_int1=D_int1*(delta_t)/(2*delta_x^2);% Calculate sigma of the first interface
sigma_int2=D_int2*(delta_t)/(2*delta_x^2);% Calculate sigma of the second interface
h=Dk(3)/Dv(1); %Calculate h

%% Build simulation components
% Set up a sub-matrix L which contains the information regarding the surface exchange as well as diffusion through the LSCF layer

% Set up the vector for the subdiagonal
sub = zeros(round(steps_x)-1,1);
sub(1:round(L(1)/delta_x)-1,1)=sigma_1;
sub(round(L(1)/delta_x),1)=sigma_int1;
sub(round(L(1)/delta_x)+1:round((L(1)+L(2))/delta_x-1),1)=sigma_2;
sub(round((L(1)+L(2))/delta_x),1)=sigma_int2;
sub(round((L(1)+L(2))/delta_x+1):end,1)=sigma_3;
sub(1:length(Dv))=sigma_v;


% Set up the vector for the superdiagonal
sup = zeros(round(steps_x)-1,1);
sup(1:round(L(1)/delta_x),1)=sigma_1;
sup(round(L(1)/delta_x),1)=sigma_int1;
sup(round(L(1)/delta_x)+1:round((L(1)+L(2))/delta_x),1)=sigma_2;
sup(round((L(1)+L(2))/delta_x),1)=sigma_int2;
sup(round((L(1)+L(2))/delta_x)+1:end,1)=sigma_3;
sup(1:length(Dv))=sigma_v;


% Set up the vector for the central diagonal
s1 = zeros(round(steps_x),1);
s1(1:round(L(1)/delta_x)-1,1)=1-2*sigma_1;
s1(round(L(1)/delta_x),1)=1-(sigma_int1+sigma_1);
s1(round(L(1)/delta_x)+1,1)=1-(sigma_int1+sigma_2);
s1(round(L(1)/delta_x)+2:round((L(1)+L(2))/delta_x)-1,1)=1-2*sigma_2;
s1(round((L(1)+L(2))/delta_x),1)=1-(sigma_2+sigma_int2);
s1(round((L(1)+L(2))/delta_x)+1,1)=1-(sigma_3+sigma_int2);
s1(round((L(1)+L(2))/delta_x)+2:end,1)=1-2*sigma_3;
s1(1)=1-(sigma_v(1)+sigma_v(2));
s1(2)=1-(sigma_v(2)+sigma_v(3));
s1(3)=1-(sigma_v(3)+sigma_v(4));
s1(4)=1-(sigma_v(4)+sigma_v(5));
s1(5)=1-(sigma_v(5)+sigma_1);


s1_n = zeros(round(steps_x),1);
s1_n(1:round(L(1)/delta_x)-1,1)=1+2*sigma_1;
s1_n(round(L(1)/delta_x),1)=1+(sigma_int1+sigma_1);
s1_n(round(L(1)/delta_x)+1,1)=1+(sigma_int1+sigma_2);
s1_n(round(L(1)/delta_x)+2:round((L(1)+L(2))/delta_x)-1,1)=1+2*sigma_2;
s1_n(round((L(1)+L(2))/delta_x),1)=1+(sigma_2+sigma_int2);
s1_n(round((L(1)+L(2))/delta_x)+1,1)=1+(sigma_3+sigma_int2);
s1_n(round((L(1)+L(2))/delta_x)+2:end,1)=1+2*sigma_3;
s1_n(1)=1+(sigma_v(1)+sigma_v(2));
s1_n(2)=1+(sigma_v(2)+sigma_v(3));
s1_n(3)=1+(sigma_v(3)+sigma_v(4));
s1_n(4)=1+(sigma_v(4)+sigma_v(5));
s1_n(5)=1+(sigma_v(5)+sigma_1);

% Create the matrix at the current timestep (Mirror boundary condition)
A = full(gallery('tridiag',sub,s1,sup));
A(1,1) = 1-2*sigma_v(1)*(1+delta_x*h); %Surface boundary condition
A(1,2) = 2*sigma_v(1); %Surface boundary condition
A(end,end) = 1+sigma_3; %Mirror boundary condition parameters
A(end,end-1) = -2*sigma_3;
A(end,end-2) = sigma_3;

% Create the matrix at the future timestep
A_n = full(gallery('tridiag',-sub,s1_n,-sup));
A_n(1,1) = 1+2*sigma_v(1)*(1+delta_x*h); %Surface boundary condition
A_n(1,2) = -2*sigma_v(1); %Surface boundary condition
A_n(end,end) = 1-sigma_3; %Mirror boundary condition parameters
A_n(end,end-1) = 2*sigma_3;
A_n(end,end-2) = -sigma_3;

% Create the matrix at the current timestep (Dirichlet boundary condition)
A_Di = full(gallery('tridiag',sub,s1,sup));
A_Di(1,1) = 1-2*sigma_v(1)*(1+delta_x*h); %Surface boundary condition
A_Di(1,2) = 2*sigma_v(1); %Surface boundary condition
A_Di(end,end) = 1; % Dirichlet boundary condition
A_Di(end,end-1) = 0;
A_Di(end,end-2) = 0;

A_n_Di = full(gallery('tridiag',-sub,-s1+2,-sup));
A_n_Di(1,1) = 1+2*sigma_v(1)*(1+delta_x*h); %Surface boundary condition
A_n_Di(1,2) = -2*sigma_v(1); %Surface boundary condition
A_n_Di(end,end) = 1; %Mirror boundary condition parameters
A_n_Di(end,end-1) = 0;
A_n_Di(end,end-2) = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Define surface exchange vector
G = zeros(round(steps_x), 1);
G(1) = 4*delta_x*sigma_v(1)*h; %surface boundary condition vector G
C = zeros(round(steps_x),1);%Initial position vector, 0 everywhere
C_Di = zeros(round(steps_x),1);% Initial position vector for the Dirichlet bc


%% Run the iterative matrix calculation

for t=1:steps_t
    C = A_n\(A*C+G);
    C_Di = A_n_Di\(A_Di*C_Di+G);
    
end

    C = (C+C_Di)/2;
    xq = 0:1000*sum(L)/(length(data)-1):1000*sum(L);
    
    C = interp1(x*1000,C,xq);
end
