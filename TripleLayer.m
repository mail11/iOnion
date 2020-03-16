%This script produces a Crank-Nicolson simulation of a stack of two layers
%of material. The main purpose is to simulate a diffusion profile of a system with a GDC
%electrolyte and an LSCF cathode as well as an interface between the two.

clear;close all

%% Define parameters
% Material properties
D_1=0.013769; % Define diffusivity of cathode material [um^2/s] (LSCF)
D_2=0.0021;  % Define diffusivity of electrolyte material [um^2/s] (GDC)
D_3=0.01; % Define diffusivity of the third layer [um^2/s] (YSZ)
k=0.0013564; % Define surface exchange coefficient of LSCF [um/s]
D_int1=0.00155; %Define interfacial resistance [s/um]
D_int2=0.01; %Define interfacial resistance at the second interface [s/um]
int_width=0;% Define the width of the interface region (0 is the default value)

% Experimental setup
Duration=.083; % Time of exchange in hours (found with Kiloran correction)
L_1 = 0.2; % Length of LSCF [um]
L_2 = 1.1; % Length of GDC [um]
L_3 = 1.1; % Length of YSZ[um]

% Simulation parameters
delta_x=0.01; % Define the spatial step [um]
delta_t=2; % Define the time step [s] (this should be a derived parameter based on max value of sigma)

% Derived properties
%D_int=2*(r/delta_x+1/D_1+1/D_2)^-1; % Define the diffusivity of the interface [um^2/s]
Length=L_1+L_2+L_3; % Total sample length [um]
steps_x=Length/delta_x+1; % Define number of spatial nodes
steps_t=Duration*3600/delta_t+1; % Define number of time steps
x=0:delta_x:Length; % Define the length of the profile in steps of delta_x
sigma_1=D_1*(delta_t)/(2*delta_x^2); % Calculate sigma of LSCF
sigma_2=D_2*(delta_t)/(2*delta_x^2); % Calculate sigma of GDC
sigma_3=D_3*(delta_t)/(2*delta_x^2); % Calculate sigma of YSZ
sigma_int1=D_int1*(delta_t)/(2*delta_x^2);% Calculate sigma of the first interface
sigma_int2=D_int2*(delta_t)/(2*delta_x^2);% Calculate sigma of the second interface
h=k/D_1; %Calculate h

%% Build simulation components
% Set up a sub-matrix L which contains the information regarding the surface exchange as well as diffusion through the LSCF layer
Layer1 = full(gallery('tridiag',round(L_1/delta_x),sigma_1,1-2*sigma_1,sigma_1)); % Bulk diffusion parameters
Layer1(1,1) = 1-2*sigma_1*(1+delta_x*h); %Surface boundary condition
Layer1(1,2) = 2*sigma_1; %Surface boundary condition

% Set up a sub-matrix E for the bulk diffusion properties of the electrolyte
Layer2 = full(gallery('tridiag',round((L_2-delta_x)/delta_x),sigma_2,1-2*sigma_2,sigma_2)); %Bulk diffusion parameters

% Set up a sub-matrix S for the bulk diffusion properties of the substrate
Layer3 = full(gallery('tridiag',round((L_3-delta_x)/delta_x)+2,sigma_3,1-2*sigma_3,sigma_3)); %Bulk diffusion parameters
Layer3(end,end) = 1+sigma_3; %Mirror boundary condition parameters
Layer3(end,end-1) = -2*sigma_3;
Layer3(end,end-2) = sigma_3;

% Set up the main matrix which combines L, E, S and interface properties into one
A = full(gallery('tridiag',round(steps_x),sigma_2,1-2*sigma_2,sigma_2)); % Create a tridiagonal basis using the electrolyte values
A(1:L_1/delta_x,1:L_1/delta_x)=Layer1; %integrate the sub-matrix L into the main matrix
A(1+(L_1/delta_x):(L_2+L_1)/delta_x-1,1+(L_1/delta_x):((L_2+L_1)/delta_x)-1)=Layer2; %integrate the sub-matrix E into the main matrix
A(round((Length-L_3)/delta_x)+1:end,round((Length-L_3)/delta_x)+1:end)=Layer3; %integrate the sub-matrix S into the main matrix

%Set up the part of the matrix at the interface with varying width
for l1=L_1/delta_x:L_1/delta_x+int_width % Creates a loop which ensures that the tridiagonal matrix at the interface is set up with varying interface width
    A(l1,l1-1) = sigma_int1;
    A(l1,l1) = 1-2*sigma_int1;
    A(l1,l1+1) = sigma_int1;
end

A(L_1/delta_x+int_width,L_1/delta_x+1+int_width) = sigma_2;
A(L_1/delta_x+1+int_width,L_1/delta_x+int_width) = sigma_int1; % Manual modification to ensure the matrix is set up correctly


%Set up the part of the matrix at the second interface with varying width
for l2=(L_1+L_2)/delta_x:(L_1+L_2)/delta_x+int_width % Creates a loop which ensures that the tridiagonal matrix at the interface is set up with varying interface width
    A(l2,l2-1) = sigma_int2;
    A(l2,l2) = 1-2*sigma_int2;
    A(l2,l2+1) = sigma_int2;
end

A((L_1+L_2)/delta_x+int_width,(L_2+L_1)/delta_x+1+int_width) = sigma_3;
A((L_1+L_2)/delta_x+1+int_width,(L_2+L_1)/delta_x+int_width) = sigma_int2; % Manual modification to ensure the matrix is set up correctly

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Setting up the transformation matrix for the next timestep A(n+1), which
%is just the same as the previous matrix but with flipped signs.

% Set up a sub-matrix L which contains the information regarding the surface exchange as well as diffusion through the LSCF layer
Layer1_n = full(gallery('tridiag',round(L_1/delta_x),-sigma_1,1+2*sigma_1,-sigma_1)); % Bulk diffusion parameters
Layer1_n(1,1) = 1+2*sigma_1*(1+delta_x*h); %Surface boundary condition
Layer1_n(1,2) = -2*sigma_1; %Surface boundary condition

% Set up a sub-matrix E for the bulk diffusion properties of the electrolyte
Layer2_n = full(gallery('tridiag',round((L_2-delta_x)/delta_x),-sigma_2,1+2*sigma_2,-sigma_2)); %Bulk diffusion parameters

% Set up a sub-matrix S for the bulk diffusion properties of the substrate
Layer3_n = full(gallery('tridiag',round((L_3-delta_x)/delta_x)+2,-sigma_3,1+2*sigma_3,-sigma_3)); %Bulk diffusion parameters
Layer3_n(end,end) = 1-sigma_3; %Mirror boundary condition parameters
Layer3_n(end,end-1) = 2*sigma_3;
Layer3_n(end,end-2) = -sigma_3;

% Set up the main matrix which combines L, E, S and interface properties into one
A_n = full(gallery('tridiag',round(steps_x),-sigma_2,1+2*sigma_2,-sigma_2)); % Create a tridiagonal basis using the electrolyte values
A_n(1:L_1/delta_x,1:L_1/delta_x)=Layer1_n; %integrate the sub-matrix L into the main matrix
A_n(1+(L_1/delta_x):(L_2+L_1)/delta_x-1,1+(L_1/delta_x):((L_2+L_1)/delta_x)-1)=Layer2_n; %integrate the sub-matrix E into the main matrix
A_n(round((Length-L_3)/delta_x)+1:end,round((Length-L_3)/delta_x)+1:end)=Layer3_n; %integrate the sub-matrix S into the main matrix

%Set up the part of the matrix at the interface with varying width
for l1=L_1/delta_x:L_1/delta_x+int_width % Creates a loop which ensures that the tridiagonal matrix at the interface is set up with varying interface width
    A_n(l1,l1-1) = -sigma_int1;
    A_n(l1,l1) = 1+2*sigma_int1;
    A_n(l1,l1+1) = -sigma_int1;
end

A_n(L_1/delta_x+int_width,L_1/delta_x+1+int_width) = -sigma_2;
A_n(L_1/delta_x+1+int_width,L_1/delta_x+int_width) = -sigma_int1; % Manual modification to ensure the matrix is set up correctly


%Set up the part of the matrix at the second interface with varying width
for l2=(L_2+L_1)/delta_x:(L_2+L_1)/delta_x+int_width % Creates a loop which ensures that the tridiagonal matrix at the interface is set up with varying interface width
    A_n(l2,l2-1) = -sigma_int2;
    A_n(l2,l2) = 1+2*sigma_int2;
    A_n(l2,l2+1) = -sigma_int2;
end

A_n((L_1+L_2)/delta_x+int_width,(L_2+L_1)/delta_x+1+int_width) = -sigma_3;
A_n((L_1+L_2)/delta_x+1+int_width,(L_2+L_1)/delta_x+int_width) = -sigma_int2; % Manual modification to ensure the matrix is set up correctly

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Define surface exchange vector
G = zeros(round(steps_x), 1);
G(1) = 4*delta_x*sigma_1*h; %surface boundary condition vector G

C = zeros(round(steps_x),1);%Initial position vector, 0 everywhere

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initialise figure
Fig=figure(...
    'Units','normalized',...
    'Position',[.3 .5 .5 .4],...
    'Color',[1 1 1],...%    'renderer','painters',...
    'WindowStyle','normal',...
    'PaperPositionMode','auto',...
    'PaperOrientation','landscape');
H=plot(x,C,'x-',[L_1,L_1],[0,1],'k:','linewidth',1.5);
set(gca,'TickLabelInterpreter','latex',...
    'LineWidth',1.2,...
    'FontSize',16);
xlim([0 Length])
ylim([0 0.5]); %inf is quite annoying to watch :-)
xlabel('Distance, $x$ / microns','interpreter','latex')
ylabel('Normalised isotopic fraction, $C''$','interpreter','latex')
linkdata on %This updates the data in the plot
set(gca,'TickLabelInterpreter','latex',...
    'LineWidth',1.2,...
    'FontSize',16);
legend('C-N Simulation','Interface location', 'Interpreter','Latex');
%% Run the iterative matrix calculation

for t=1:steps_t
    C = A_n\(A*C+G);
    
    % Update plot
    refreshdata
    drawnow
end

hold off

% Optional code to produce noisy data that can be used to test fitting
% functions

% snr = 50;
% out = awgn(C,snr);varargin;
%
% plot(out)
%
% writematrix(out,'GenNoiseInt.txt')
%