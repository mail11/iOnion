%This script produces a Crank-Nicolson simulation of a stack of two layers
%of material. The main purpose is to simulate a diffusion profile of a system with a GDC
%electrolyte and an LSCF cathode as well as an interface between the two.

clear;close all

%% Define parameters
% Material properties
D_LSCF=.29; % Define diffusivity of cathode material [um^2/s] (LSCF)
D_GDC=.40;  % Define diffusivity of electrolyte material [um^2/s] (GDC)
k=0.00772; % Define surface exchange coefficient of LSCF [um/s]
D_int=.30; %Define interfacial resistance [s/um]
int_width=0;% Define the width of the interface region (0 is the default value)

% Experimental setup
Duration=.05; % Time of exchange in hours (found with Kiloran correction)
L_1 = 0.2; % Length of LSCF [um]
L_2 = 1; % Length of GDC [um]

% Simulation parameters
delta_x=0.01; % Define the spatial step [um]
delta_t=2; % Define the time step [s] (this should be a derived parameter based on max value of sigma)

% Derived properties
%D_int=2*(r/delta_x+1/D_LSCF+1/D_GDC)^-1; % Define the diffusivity of the interface [um^2/s]
Length=L_1+L_2; % Total sample length [um]
steps_x=Length/delta_x+1; % Define number of spatial nodes
steps_t=Duration*3600/delta_t+1; % Define number of time steps
x=0:delta_x:Length; % Define the length of the profile in steps of delta_x
sigmaLSCF=D_LSCF*(delta_t)/(2*delta_x^2); % Calculate sigma of LSCF
sigmaGDC=D_GDC*(delta_t)/(2*delta_x^2); % Calculate sigma of GDC
sigma_int=D_int*(delta_t)/(2*delta_x^2);% Calculate sigma of the first interface
h=k/D_LSCF; %Calculate h

%% Build simulation components
% Set up a sub-matrix L which contains the information regarding the surface exchange as well as diffusion through the LSCF layer
Layer1 = full(gallery('tridiag',round(L_1/delta_x),sigmaLSCF,1-2*sigmaLSCF,sigmaLSCF)); % Bulk diffusion parameters
Layer1(1,1) = 1-2*sigmaLSCF*(1+delta_x*h); %Surface boundary condition
Layer1(1,2) = 2*sigmaLSCF; %Surface boundary condition

% Set up a sub-matrix E for the bulk diffusion properties of the electrolyte
Layer2 = full(gallery('tridiag',round((L_2-delta_x)/delta_x)+2,sigmaGDC,1-2*sigmaGDC,sigmaGDC)); %Bulk diffusion parameters
Layer2(end,end) = 1+sigmaGDC; %Mirror boundary condition parameters
Layer2(end,end-1) = -2*sigmaGDC;
Layer2(end,end-2) = sigmaGDC;

% Set up the main matrix which combines L, E and interface properties into one
A = full(gallery('tridiag',round(steps_x),sigmaGDC,1-2*sigmaGDC,sigmaGDC)); % Create a tridiagonal basis using the electrolyte values
A(1:L_1/delta_x,1:L_1/delta_x)=Layer1; %integrate the sub-matrix L into the main matrix
A(round((Length-L_2)/delta_x)+1:end,round((Length-L_2)/delta_x)+1:end)=Layer2; %integrate the sub-matrix E into the main matrix

%Set up the part of the matrix at the interface with varying width
for l=L_1/delta_x:L_1/delta_x+int_width % Creates a loop which ensures that the tridiagonal matrix at the interface is set up with varying interface width
    A(l,l-1) = sigma_int;
    A(l,l) = 1-2*sigma_int;
    A(l,l+1) = sigma_int;
end

A(L_1/delta_x+int_width,L_1/delta_x+1+int_width) = sigmaGDC;
A(L_1/delta_x+1+int_width,L_1/delta_x+int_width) = sigma_int; % Manual modification to ensure the matrix is set up correctly

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Setting up the transformation matrix for the next timestep A(n+1), which
%is just the same as the previous matrix but with flipped signs.

Layer1_n = full(gallery('tridiag',round(L_1/delta_x),-sigmaLSCF,1+2*sigmaLSCF,-sigmaLSCF));
Layer1_n(1,1) = 1+2*sigmaLSCF*(1+delta_x*h);
Layer1_n(1,2) = -2*sigmaLSCF;

Layer2_n = full(gallery('tridiag',round((L_2-delta_x)/delta_x)+2,-sigmaGDC,1+2*sigmaGDC,-sigmaGDC));
Layer2_n(end,end) = 1-sigmaGDC;
Layer2_n(end,end-1) = 2*sigmaGDC;
Layer2_n(end,end-2) = -sigmaGDC;

A_n = full(gallery('tridiag',round(steps_x),-sigmaGDC,1+2*sigmaGDC,-sigmaGDC));
A_n(1:L_1/delta_x,1:L_1/delta_x)=Layer1_n;
A_n(round((Length-L_2)/delta_x)+1:end,round((Length-L_2)/delta_x)+1:end)=Layer2_n;

for l=L_1/delta_x:(L_1+int_width)/delta_x
    A_n(l,l-1) = -sigma_int;
    A_n(l,l) = 1+2*sigma_int;
    A_n(l,l+1) = -sigma_int;
end

A_n(L_1/delta_x+int_width,L_1/delta_x+1+int_width) = -sigmaGDC;
A_n(L_1/delta_x+1+int_width,L_1/delta_x+int_width) = -sigma_int;

% Define surface exchange vector
G = zeros(round(steps_x), 1);
G(1) = 4*delta_x*sigmaLSCF*h; %surface boundary condition vector G

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
h=plot(x,C,'x-',[L_1,L_1],[0,1],'k:','linewidth',1.5);
set(gca,'TickLabelInterpreter','latex',...
    'LineWidth',1.2,...
    'FontSize',16);
xlim([0 Length])
ylim([0 1]); %inf is quite annoying to watch :-)
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