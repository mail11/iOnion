function C = CNInt(Dk)

%% Define parameters

% Set up the global variables

global data;
global L;
global delta_x;
global steps_x;
global Duration;
global Layers;

int_width=0;% Define the width of the interface region (0 is the default value)

% Simulation parameters
delta_t=2; % Define the time step [s] (this should be a derived parameter based on max value of sigma)

if Layers == 1
    
    %% One layer
    
    Length=sum(L); % Total sample length [um]
    delta_x=Length/(length(data)-1); % Define the spatial step [um]
    steps_x=Length/delta_x+1; % Define number of spatial nodes
    steps_t=Duration*3600/delta_t+1; % Define number of time steps
    x=0:delta_x:Length; % Define the length of the profile in steps of delta_x
    sigma_1=(10^(Dk(1)))*(delta_t)/(2*delta_x^2); % Calculate sigma of LSCF
    h=(10^(Dk(2)))/(10^(Dk(1))); %Calculate h
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Build simulation components
    sub = zeros(round(steps_x)-1,1);
    sub(1:end,1)=sigma_1;
    
    sup = zeros(round(steps_x)-1,1);
    sup(1:end,1)=sigma_1;
    
    s1 = zeros(round(steps_x),1);
    s1(1:end,1)=1-2*sigma_1;
    
    A = full(gallery('tridiag',sub,s1,sup));
    A(1,1) = 1-2*sigma_1*(1+delta_x*h); %Surface boundary condition
    A(1,2) = 2*sigma_1; %Surface boundary condition
    A(end,end) = 1+sigma_1; %Mirror boundary condition parameters
    A(end,end-1) = -2*sigma_1;
    A(end,end-2) = sigma_1;
    
    A_n = full(gallery('tridiag',-sub,-s1+2,-sup));
    A_n(1,1) = 1+2*sigma_1*(1+delta_x*h); %Surface boundary condition
    A_n(1,2) = -2*sigma_1; %Surface boundary condition
    A_n(end,end) = 1-sigma_1; %Mirror boundary condition parameters
    A_n(end,end-1) = 2*sigma_1;
    A_n(end,end-2) = -sigma_1;
    
    A_Di = full(gallery('tridiag',sub,s1,sup));
    A_Di(1,1) = 1-2*sigma_1*(1+delta_x*h); %Surface boundary condition
    A_Di(1,2) = 2*sigma_1; %Surface boundary condition
    A_Di(end,end) = 1; % Dirichlet boundary condition
    A_Di(end,end-1) = 0;
    A_Di(end,end-2) = 0;
    
    A_n_Di = full(gallery('tridiag',-sub,-s1+2,-sup));
    A_n_Di(1,1) = 1+2*sigma_1*(1+delta_x*h); %Surface boundary condition
    A_n_Di(1,2) = -2*sigma_1; %Surface boundary condition
    A_n_Di(end,end) = 1; %Mirror boundary condition parameters
    A_n_Di(end,end-1) = 0;
    A_n_Di(end,end-2) = 0;
    
elseif Layers == 2 % If there are two layers, fminsearch will call upon the following code which is written for a matrix containing two layers
    %% Two layers
    
    % Derived properties
    %D_int=2*(r/delta_x+1/D_LSCF+1/D_GDC)^-1; % Define the diffusivity of the interface [um^2/s]
    Length=sum(L); % Total sample length [um]
    delta_x=Length/(length(data)-1); % Define the spatial step [um]
    steps_x=Length/delta_x+1; % Define number of spatial nodes
    steps_t=Duration*3600/delta_t+1; % Define number of time steps
    x=0:delta_x:Length; % Define the length of the profile in steps of delta_x
    sigma_1=(10^(Dk(1)))*(delta_t)/(2*delta_x^2); % Calculate sigma of LSCF
    sigma_2=(10^(Dk(2)))*(delta_t)/(2*delta_x^2); % Calculate sigma of GDC
    sigma_int=(10^(Dk(3)))*(delta_t)/(2*delta_x^2);% Calculate sigma of the first interface
    h=(10^(Dk(4)))/(10^(Dk(1))); %Calculate h
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Build simulation components
    % Set up the vector for the subdiagonal
    sub = zeros(round(steps_x)-1,1);
    sub(1:round(L(1)/delta_x)-1,1)=sigma_1;
    sub(round(L(1)/delta_x):round(L(1)/delta_x)+1,1)=sigma_int;
    sub(round(L(1)/delta_x)+1:end,1)=sigma_2;
    
    % Set up the vector for the superdiagonal
    sup = zeros(round(steps_x)-1,1);
    sup(1:round(L(1)/delta_x)-1,1)=sigma_1;
    sup(round(L(1)/delta_x)-0:round(L(1)/delta_x),1)=sigma_int;
    sup(round(L(1)/delta_x)+1:end,1)=sigma_2;
    
    % Set up the vector for the central diagonal
    s1 = zeros(round(steps_x),1);
    s1(1:round(L(1)/delta_x)-1,1)=1-2*sigma_1;
    s1(round(L(1)/delta_x),1)=1-(sigma_int+sigma_1);
    s1(round(L(1)/delta_x)+1,1)=1-(sigma_int+sigma_2);
    s1(round(L(1)/delta_x)+2:end,1)=1-2*sigma_2;
    
    s1_n = zeros(round(steps_x),1);
    s1_n(1:round(L(1)/delta_x)-1,1)=1+2*sigma_1;
    s1_n(round(L(1)/delta_x),1)=1+(sigma_int+sigma_1);
    s1_n(round(L(1)/delta_x)+1,1)=1+(sigma_int+sigma_2);
    s1_n(round(L(1)/delta_x)+2:end,1)=1+2*sigma_2;
    
    % Create the matrix at the current timestep (Mirror boundary condition)
    A = full(gallery('tridiag',sub,s1,sup));
    A(1,1) = 1-2*sigma_1*(1+delta_x*h); %Surface boundary condition
    A(1,2) = 2*sigma_1; %Surface boundary condition
    A(end,end) = 1+sigma_2; %Mirror boundary condition parameters
    A(end,end-1) = -2*sigma_2;
    A(end,end-2) = sigma_2;
    
    % Create the matrix at the future timestep
    A_n = full(gallery('tridiag',-sub,s1_n,-sup));
    A_n(1,1) = 1+2*sigma_1*(1+delta_x*h); %Surface boundary condition
    A_n(1,2) = -2*sigma_1; %Surface boundary condition
    A_n(end,end) = 1-sigma_2; %Mirror boundary condition parameters
    A_n(end,end-1) = 2*sigma_2;
    A_n(end,end-2) = -sigma_2;
    
    % Create the matrix at the current timestep (Dirichlet boundary condition)
    A_Di = full(gallery('tridiag',sub,s1,sup));
    A_Di(1,1) = 1-2*sigma_1*(1+delta_x*h); %Surface boundary condition
    A_Di(1,2) = 2*sigma_1; %Surface boundary condition
    A_Di(end,end) = 1; % Dirichlet boundary condition
    A_Di(end,end-1) = 0;
    A_Di(end,end-2) = 0;
    
    A_n_Di = full(gallery('tridiag',-sub,-s1+2,-sup));
    A_n_Di(1,1) = 1+2*sigma_1*(1+delta_x*h); %Surface boundary condition
    A_n_Di(1,2) = -2*sigma_1; %Surface boundary condition
    A_n_Di(end,end) = 1; %Mirror boundary condition parameters
    A_n_Di(end,end-1) = 0;
    A_n_Di(end,end-2) = 0;
    
elseif Layers == 3 % Set up a matric calculation with three layers
    
    %% Three layers
    
    % Derived properties
    %D_int=2*(r/delta_x+1/(10^(Dk(1)))+1/(10^(Dk(2))))^-1; % Define the diffusivity of the interface [um^2/s]
    Length=sum(L); % Total sample length [um]
    delta_x=Length/(length(data)-1);
    steps_x=Length/delta_x+1; % Define number of spatial nodes
    steps_t=Duration*3600/delta_t+1; % Define number of time steps
    x=0:delta_x:Length; % Define the length of the profile in steps of delta_x
    sigma_1=(10^(Dk(1)))*(delta_t)/(2*delta_x^2); % Calculate sigma of LSCF
    sigma_2=(10^(Dk(2)))*(delta_t)/(2*delta_x^2); % Calculate sigma of GDC
    sigma_3=(10^(Dk(3)))*(delta_t)/(2*delta_x^2); % Calculate sigma of YSZ
    sigma_int1=(10^(Dk(4)))*(delta_t)/(2*delta_x^2);% Calculate sigma of the first interface
    sigma_int2=(10^(Dk(5)))*(delta_t)/(2*delta_x^2);% Calculate sigma of the second interface
    h=(10^(Dk(6)))/(10^(Dk(1))); %Calculate h
    
    
    %% Build simulation components
    % Set up a sub-matrix L which contains the information regarding the surface exchange as well as diffusion through the LSCF layer
    % Set up the vector for the subdiagonal
    % Set up the vector for the subdiagonal
    sub = zeros(round(steps_x)-1,1);
    sub(1:round(L(1)/delta_x)-1,1)=sigma_1;
    sub(round(L(1)/delta_x),1)=sigma_int1;
    sub(round(L(1)/delta_x)+1:round((L(1)+L(2))/delta_x-1),1)=sigma_2;
    sub(round((L(1)+L(2))/delta_x),1)=sigma_int2;
    sub(round((L(1)+L(2))/delta_x+1):end,1)=sigma_3;
    
    % Set up the vector for the superdiagonal
    sup = zeros(round(steps_x)-1,1);
    sup(1:round(L(1)/delta_x),1)=sigma_1;
    sup(round(L(1)/delta_x),1)=sigma_int1;
    sup(round(L(1)/delta_x)+1:round((L(1)+L(2))/delta_x),1)=sigma_2;
    sup(round((L(1)+L(2))/delta_x),1)=sigma_int2;
    sup(round((L(1)+L(2))/delta_x)+1:end,1)=sigma_3;
    
    % Set up the vector for the central diagonal
    s1 = zeros(round(steps_x),1);
    s1(1:round(L(1)/delta_x)-1,1)=1-2*sigma_1;
    s1(round(L(1)/delta_x),1)=1-(sigma_int1+sigma_1);
    s1(round(L(1)/delta_x)+1,1)=1-(sigma_int1+sigma_2);
    s1(round(L(1)/delta_x)+2:round((L(1)+L(2))/delta_x)-1,1)=1-2*sigma_2;
    s1(round((L(1)+L(2))/delta_x),1)=1-(sigma_2+sigma_int2);
    s1(round((L(1)+L(2))/delta_x)+1,1)=1-(sigma_3+sigma_int2);
    s1(round((L(1)+L(2))/delta_x)+2:end,1)=1-2*sigma_3;
    
    
    s1_n = zeros(round(steps_x),1);
    s1_n(1:round(L(1)/delta_x)-1,1)=1+2*sigma_1;
    s1_n(round(L(1)/delta_x),1)=1+(sigma_int1+sigma_1);
    s1_n(round(L(1)/delta_x)+1,1)=1+(sigma_int1+sigma_2);
    s1_n(round(L(1)/delta_x)+2:round((L(1)+L(2))/delta_x)-1,1)=1+2*sigma_2;
    s1_n(round((L(1)+L(2))/delta_x),1)=1+(sigma_2+sigma_int2);
    s1_n(round((L(1)+L(2))/delta_x)+1,1)=1+(sigma_3+sigma_int2);
    s1_n(round((L(1)+L(2))/delta_x)+2:end,1)=1+2*sigma_3;
    
    % Create the matrix at the current timestep (Mirror boundary condition)
    A = full(gallery('tridiag',sub,s1,sup));
    A(1,1) = 1-2*sigma_1*(1+delta_x*h); %Surface boundary condition
    A(1,2) = 2*sigma_1; %Surface boundary condition
    A(end,end) = 1+sigma_3; %Mirror boundary condition parameters
    A(end,end-1) = -2*sigma_3;
    A(end,end-2) = sigma_3;
    
    % Create the matrix at the future timestep
    A_n = full(gallery('tridiag',-sub,s1_n,-sup));
    A_n(1,1) = 1+2*sigma_1*(1+delta_x*h); %Surface boundary condition
    A_n(1,2) = -2*sigma_1; %Surface boundary condition
    A_n(end,end) = 1-sigma_3; %Mirror boundary condition parameters
    A_n(end,end-1) = 2*sigma_3;
    A_n(end,end-2) = -sigma_3;
    
    %Create the matrix at the current timestep (Dirichlet boundary condition)
    A_Di = full(gallery('tridiag',sub,s1,sup));
    A_Di(1,1) = 1-2*sigma_1*(1+delta_x*h); %Surface boundary condition
    A_Di(1,2) = 2*sigma_1; %Surface boundary condition
    A_Di(end,end) = 1; % Dirichlet boundary condition
    A_Di(end,end-1) = 0;
    A_Di(end,end-2) = 0;
    
    A_n_Di = full(gallery('tridiag',-sub,-s1+2,-sup));
    A_n_Di(1,1) = 1+2*sigma_1*(1+delta_x*h); %Surface boundary condition
    A_n_Di(1,2) = -2*sigma_1; %Surface boundary condition
    A_n_Di(end,end) = 1; %Mirror boundary condition parameters
    A_n_Di(end,end-1) = 0;
    A_n_Di(end,end-2) = 0;
end

% Define surface exchange vector
G = zeros(round(steps_x), 1);
G(1) = 4*delta_x*sigma_1*h; %surface boundary condition vector G

C = zeros(round(steps_x),1);%Initial position vector, 0 everywhere
C_Di = zeros(round(steps_x),1);

%% Run the iterative matrix calculation

for t=1:steps_t
    
    C = A_n\(A*C+G);
    C_Di = A_n_Di\(A_Di*C_Di+G); % Calculate Dirichlet boundary condition
    %if desired
    
end

C = (C+C_Di)/2;
end