function C = CNInt(Dk)

%% Define parameters

% Set up the global variables

global data;
global L_1;
global L_2;
global L_3;
global steps_x;
global answer;
global Duration;

int_width=0;% Define the width of the interface region (0 is the default value)

% Simulation parameters
delta_t=2; % Define the time step [s] (this should be a derived parameter based on max value of sigma)

if str2double(answer(2)) == 1
    
    %% One layer
    
    Length=L_1; % Total sample length [um]
    delta_x=Length/(length(data)-1); % Define the spatial step [um]
    steps_x=Length/delta_x+1; % Define number of spatial nodes
    steps_t=Duration*3600/delta_t+1; % Define number of time steps
    x=0:delta_x:Length; % Define the length of the profile in steps of delta_x
    sigma_1=(10^(Dk(1)))*(delta_t)/(2*delta_x^2); % Calculate sigma of LSCF
    h=(10^(Dk(2)))/(10^(Dk(1))); %Calculate h
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Build simulation components
    % Set up a sub-matrix L which contains the information regarding the surface exchange as well as diffusion through the LSCF layer
    A = full(gallery('tridiag',round(steps_x),sigma_1,1-2*sigma_1,sigma_1)); % Create a tridiagonal basis using the electrolyte values
    A(1,1) = 1-2*sigma_1*(1+delta_x*h); %Surface boundary condition
    A(1,2) = 2*sigma_1; %Surface boundary condition
    A(end,end) = 1+sigma_1; %Mirror boundary condition parameters
    A(end,end-1) = -2*sigma_1;
    A(end,end-2) = sigma_1;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %Setting up the transformation matrix for the next timestep A(n+1), which
    %is just the same as the previous matrix but with flipped signs.
    
    A_n = full(gallery('tridiag',round(steps_x),-sigma_1,1+2*sigma_1,-sigma_1));
    A_n(1,1) = 1+2*sigma_1*(1+delta_x*h);
    A_n(1,2) = -2*sigma_1;
    A_n(end,end) = 1-sigma_1;
    A_n(end,end-1) = 2*sigma_1;
    A_n(end,end-2) = -sigma_1;
    
elseif str2double(answer(2)) == 2 % If there are two layers, fminsearch will call upon the following code which is written for a matrix containing two layers
    %% Two layers
    
    % Derived properties
    %D_int=2*(r/delta_x+1/D_LSCF+1/D_GDC)^-1; % Define the diffusivity of the interface [um^2/s]
    Length=L_1+L_2; % Total sample length [um]
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
    % Set up a sub-matrix L which contains the information regarding the surface exchange as well as diffusion through the LSCF layer
    Layer1 = full(gallery('tridiag',round(L_1/delta_x),sigma_1,1-2*sigma_1,sigma_1)); % Bulk diffusion parameters
    Layer1(1,1) = 1-2*sigma_1*(1+delta_x*h); %Surface boundary condition
    Layer1(1,2) = 2*sigma_1; %Surface boundary condition
    
    % Set up a sub-matrix E for the bulk diffusion properties of the electrolyte
    Layer2 = full(gallery('tridiag',round((L_2-delta_x)/delta_x)+2,sigma_2,1-2*sigma_2,sigma_2)); %Bulk diffusion parameters
    Layer2(end,end) = 1+sigma_2; %Mirror boundary condition parameters
    Layer2(end,end-1) = -2*sigma_2;
    Layer2(end,end-2) = sigma_2;
    
    % Set up the main matrix which combines L, E and interface properties into one
    A = full(gallery('tridiag',round(steps_x),sigma_2,1-2*sigma_2,sigma_2)); % Create a tridiagonal basis using the electrolyte values
    A(1:round(L_1/delta_x),1:round(L_1/delta_x))=Layer1; %integrate the sub-matrix L into the main matrix
    A(round((Length-L_2)/delta_x)+1:end,round((Length-L_2)/delta_x)+1:end)=Layer2; %integrate the sub-matrix E into the main matrix
    
    %Set up the part of the matrix at the interface with varying width
    for l=round(L_1/delta_x):round(L_1/delta_x)+int_width % Creates a loop which ensures that the tridiagonal matrix at the interface is set up with varying interface width
        A(l,l-1) = sigma_int;
        A(l,l) = 1-2*sigma_int;
        A(l,l+1) = sigma_int;
    end
    
    A(round(L_1/delta_x)+int_width,round(L_1/delta_x)+1+int_width) = sigma_2;
    A(round(L_1/delta_x)+1+int_width,round(L_1/delta_x)+int_width) = sigma_int; % Manual modification to ensure the matrix is set up correctly
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %Setting up the transformation matrix for the next timestep A(n+1), which
    %is just the same as the previous matrix but with flipped signs.
    
    Layer1_n = full(gallery('tridiag',round(L_1/delta_x),-sigma_1,1+2*sigma_1,-sigma_1));
    Layer1_n(1,1) = 1+2*sigma_1*(1+delta_x*h);
    Layer1_n(1,2) = -2*sigma_1;
    
    Layer2_n = full(gallery('tridiag',round((L_2-delta_x)/delta_x)+2,-sigma_2,1+2*sigma_2,-sigma_2));
    Layer2_n(end,end) = 1-sigma_2;
    Layer2_n(end,end-1) = 2*sigma_2;
    Layer2_n(end,end-2) = -sigma_2;
    
    A_n = full(gallery('tridiag',round(steps_x),-sigma_2,1+2*sigma_2,-sigma_2));
    A_n(1:round(L_1/delta_x),1:round(L_1/delta_x))=Layer1_n;
    A_n(round((Length-L_2)/delta_x)+1:end,round((Length-L_2)/delta_x)+1:end)=Layer2_n;
    
    for l=round(L_1/delta_x):round((L_1+int_width)/delta_x)
        A_n(l,l-1) = -sigma_int;
        A_n(l,l) = 1+2*sigma_int;
        A_n(l,l+1) = -sigma_int;
    end
    
    A_n(round(L_1/delta_x)+int_width,round(L_1/delta_x)+1+int_width) = -sigma_2;
    A_n(round(L_1/delta_x)+1+int_width,round(L_1/delta_x)+int_width) = -sigma_int;
    
elseif str2double(answer(2)) == 3 % Set up a matric calculation with three layers
    
    %% Three layers
    
    % Derived properties
    %D_int=2*(r/delta_x+1/(10^(Dk(1)))+1/(10^(Dk(2))))^-1; % Define the diffusivity of the interface [um^2/s]
    Length=L_1+L_2+L_3; % Total sample length [um]
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
    Layer1 = full(gallery('tridiag',round(L_1/delta_x),sigma_1,1-2*sigma_1,sigma_1)); % Bulk diffusion parameters
    Layer1(1,1) = 1-2*sigma_1*(1+delta_x*h); %Surface boundary condition
    Layer1(1,2) = 2*sigma_1; %Surface boundary condition
    
    % Set up a sub-matrix E for the bulk diffusion properties of the electrolyte
    Layer2 = full(gallery('tridiag',round((L_2-delta_x)/delta_x)+1,sigma_2,1-2*sigma_2,sigma_2)); %Bulk diffusion parameters
    
    % Set up a sub-matrix S for the bulk diffusion properties of the substrate
    Layer3 = full(gallery('tridiag',round((L_3-delta_x)/delta_x)+2,sigma_3,1-2*sigma_3,sigma_3)); %Bulk diffusion parameters
    Layer3(end,end) = 1+sigma_3; %Mirror boundary condition parameters
    Layer3(end,end-1) = -2*sigma_3;
    Layer3(end,end-2) = sigma_3;
    
    % Set up the main matrix which combines L, E, S and interface properties into one
    A = full(gallery('tridiag',round(steps_x),sigma_2,1-2*sigma_2,sigma_2)); % Create a tridiagonal basis using the electrolyte values
    A(1:round(L_1/delta_x),1:round(L_1/delta_x))=Layer1; %integrate the sub-matrix L into the main matrix
    A(1+round(L_1/delta_x):round((L_2+L_1)/delta_x),1+round(L_1/delta_x):round((L_2+L_1)/delta_x))=Layer2; %integrate the sub-matrix E into the main matrix
    A(round((Length-L_3)/delta_x)+1:end,round((Length-L_3)/delta_x)+1:end)=Layer3; %integrate the sub-matrix S into the main matrix
    
    %Set up the part of the matrix at the interface with varying width
    for l1=round(L_1/delta_x):round(L_1/delta_x)+int_width % Creates a loop which ensures that the tridiagonal matrix at the interface is set up with varying interface width
        A(l1,l1-1) = sigma_int1;
        A(l1,l1) = 1-2*sigma_int1;
        A(l1,l1+1) = sigma_int1;
    end
    
    A(round(L_1/delta_x)+int_width,round(L_1/delta_x)+1+int_width) = sigma_2;
    A(round(L_1/delta_x)+1+int_width,round(L_1/delta_x)+int_width) = sigma_int1; % Manual modification to ensure the matrix is set up correctly
    
    
    %Set up the part of the matrix at the second interface with varying width
    for l2=round((L_1+L_2)/delta_x):round((L_1+L_2)/delta_x)+int_width % Creates a loop which ensures that the tridiagonal matrix at the interface is set up with varying interface width
        A(l2,l2-1) = sigma_int2;
        A(l2,l2) = 1-2*sigma_int2;
        A(l2,l2+1) = sigma_int2;
    end
    
    A(round((L_1+L_2)/delta_x)+int_width,round((L_2+L_1)/delta_x)+1+int_width) = sigma_3;
    A(round((L_1+L_2)/delta_x)+1+int_width,round((L_2+L_1)/delta_x)+int_width) = sigma_int2; % Manual modification to ensure the matrix is set up correctly
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %Setting up the transformation matrix for the next timestep A(n+1), which
    %is just the same as the previous matrix but with flipped signs.
    
    % Set up a sub-matrix L which contains the information regarding the surface exchange as well as diffusion through the LSCF layer
    Layer1_n = full(gallery('tridiag',round(L_1/delta_x),-sigma_1,1+2*sigma_1,-sigma_1)); % Bulk diffusion parameters
    Layer1_n(1,1) = 1+2*sigma_1*(1+delta_x*h); %Surface boundary condition
    Layer1_n(1,2) = -2*sigma_1; %Surface boundary condition
    
    % Set up a sub-matrix E for the bulk diffusion properties of the electrolyte
    Layer2_n = full(gallery('tridiag',round((L_2-delta_x)/delta_x+1),-sigma_2,1+2*sigma_2,-sigma_2)); %Bulk diffusion parameters
    
    % Set up a sub-matrix S for the bulk diffusion properties of the substrate
    Layer3_n = full(gallery('tridiag',round((L_3-delta_x)/delta_x)+2,-sigma_3,1+2*sigma_3,-sigma_3)); %Bulk diffusion parameters
    Layer3_n(end,end) = 1-sigma_3; %Mirror boundary condition parameters
    Layer3_n(end,end-1) = 2*sigma_3;
    Layer3_n(end,end-2) = -sigma_3;
    
    % Set up the main matrix which combines L, E, S and interface properties into one
    A_n = full(gallery('tridiag',round(steps_x),-sigma_2,1+2*sigma_2,-sigma_2)); % Create a tridiagonal basis using the electrolyte values
    A_n(1:round(L_1/delta_x),1:round(L_1/delta_x))=Layer1_n; %integrate the sub-matrix L into the main matrix
    A_n(1+round(L_1/delta_x):round((L_2+L_1)/delta_x),1+round(L_1/delta_x):(round((L_2+L_1)/delta_x)))=Layer2_n; %integrate the sub-matrix E into the main matrix
    A_n(round((Length-L_3)/delta_x)+1:end,round((Length-L_3)/delta_x)+1:end)=Layer3_n; %integrate the sub-matrix S into the main matrix
    
    %Set up the part of the matrix at the interface with varying width
    for l1=round(L_1/delta_x):round(L_1/delta_x)+int_width % Creates a loop which ensures that the tridiagonal matrix at the interface is set up with varying interface width
        A_n(l1,l1-1) = -sigma_int1;
        A_n(l1,l1) = 1+2*sigma_int1;
        A_n(l1,l1+1) = -sigma_int1;
    end
    
    A_n(round(L_1/delta_x)+int_width,round(L_1/delta_x)+1+int_width) = -sigma_2;
    A_n(round(L_1/delta_x)+1+int_width,round(L_1/delta_x)+int_width) = -sigma_int1; % Manual modification to ensure the matrix is set up correctly
    
    
    %Set up the part of the matrix at the second interface with varying width
    for l2=round((L_2+L_1)/delta_x):round((L_2+L_1)/delta_x)+int_width % Creates a loop which ensures that the tridiagonal matrix at the interface is set up with varying interface width
        A_n(l2,l2-1) = -sigma_int2;
        A_n(l2,l2) = 1+2*sigma_int2;
        A_n(l2,l2+1) = -sigma_int2;
    end
    
    A_n(round((L_1+L_2)/delta_x)+int_width,round((L_2+L_1)/delta_x)+1+int_width) = -sigma_3;
    A_n(round((L_1+L_2)/delta_x)+1+int_width,round((L_2+L_1)/delta_x)+int_width) = -sigma_int2; % Manual modification to ensure the matrix is set up correctly
    
    
end

% Define surface exchange vector
G = zeros(round(steps_x), 1);
G(1) = 4*delta_x*sigma_1*h; %surface boundary condition vector G

C = zeros(round(steps_x),1);%Initial position vector, 0 everywhere

%% Run the iterative matrix calculation

for t=1:steps_t
    
    C = A_n\(A*C+G);
    
end


end