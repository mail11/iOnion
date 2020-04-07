%This script produces a Crank-Nicolson simulation of a stack of two layers
%of material. The main purpose is to simulate a diffusion profile of a system with a GDC
%electrolyte and an LSCF cathode as well as an interface between the two.

clear;close all

%% Define parameters
% Material properties
D_1=.00026; % Define diffusivity of cathode material [um^2/s] (LSCF)
D_2=0.00026;  % Define diffusivity of electrolyte material [um^2/s] (GDC)
k=0.000292; % Define surface exchange coefficient of LSCF [um/s]
D_int=.00026; %Define interfacial resistance [s/um]
int_width=0;% Define the width of the interface region (0 is the default value)

% Experimental setup
Duration=.083; % Time of exchange in hours (found with Kiloran correction)
L = [0.2;0.8]; % Vector of the layer lengths, first element is the first layer

% Simulation parameters
delta_x=0.01; % Define the spatial step [um]
delta_t=2; % Define the time step [s] (this should be a derived parameter based on max value of sigma)

% Derived properties
%D_int=2*(r/delta_x+1/D_1+1/D_2)^-1; % Define the diffusivity of the interface [um^2/s]
Length=L(1)+L(2); % Total sample length [um]
steps_x=Length/delta_x+1; % Define number of spatial nodes
steps_t=Duration*3600/delta_t+1; % Define number of time steps
x=0:delta_x:Length; % Define the length of the profile in steps of delta_x
sigma_1=D_1*(delta_t)/(2*delta_x^2); % Calculate sigma of LSCF
sigma_2=D_2*(delta_t)/(2*delta_x^2); % Calculate sigma of GDC
sigma_int=D_int*(delta_t)/(2*delta_x^2);% Calculate sigma of the first interface
h=k/D_1; %Calculate h

%% Build simulation components
% Set up a sub-matrix L which contains the information regarding the surface exchange as well as diffusion through the LSCF layer

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Define surface exchange vector
G = zeros(round(steps_x), 1);
G(1) = 4*delta_x*sigma_1*h; %surface boundary condition vector G
C = zeros(round(steps_x),1);%Initial position vector, 0 everywhere
C_Di = zeros(round(steps_x),1);% Initial position vector for the Dirichlet bc


%% Run the iterative matrix calculation

for t=1:steps_t
    C = A_n\(A*C+G);
    C_Di = A_n_Di\(A_Di*C_Di+G);
    plot(x,C_Di,x,C,x,(C+C_Di)/2); % Plot both boundary conditions as well as their mean
    xlim([0 Length])
    ylim([-0 inf]);
    xlabel('Distance, x / microns')
    ylabel('Concentration, C')
    drawnow
end

hold off

%Optional code to produce noisy data that can be used to test fitting
%functions

% snr = 50;
% out = awgn(C,snr);varargin;
% 
% plot(out)
% 
% writematrix(out,'GenNoiseInt.txt')



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initialise figure
% Fig=figure(...
%     'Units','normalized',...
%     'Position',[.3 .5 .5 .4],...
%     'Color',[1 1 1],...%    'renderer','painters',...
%     'WindowStyle','normal',...
%     'PaperPositionMode','auto',...
%     'PaperOrientation','landscape');
% h=plot(x,C,x,C_Di,x,C_t,'x-',[L(1),L(1)],[0,1],'k:','linewidth',1.5);
% set(gca,'TickLabelInterpreter','latex',...
%     'LineWidth',1.2,...
%     'FontSize',16);
% xlim([0 Length])
% ylim([0 1]); %inf is quite annoying to watch :-)
% xlabel('Distance, $x$ / microns','interpreter','latex')
% ylabel('Normalised isotopic fraction, $C''$','interpreter','latex')
% linkdata on %This updates the data in the plot
% set(gca,'TickLabelInterpreter','latex',...
%     'LineWidth',1.2,...
%     'FontSize',16);
% legend('C-N Simulation','Interface location', 'Interpreter','Latex');
