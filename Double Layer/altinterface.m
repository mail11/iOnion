%This script produces a Crank-Nicolson simulation of a stack of two layers
%of material. The main purpose is to simulate a diffusion profile of a system with a GDC
%electrolyte and an LSCF cathode as well as an interface between the two.

close all
DLSCF=30; % Define diffusivity of cathode material [um^2/s] (LSCF)
DGDC=200; % Define diffusivity of electrolyte material [um^2/s] (GDC)
D_int=26.087; % Define the diffusivity of the interface [um^2/s]
k=0.0178; % Define surface exchange coefficient of LSCF [um/s]
Duration=1; % Time of exchange in hours
l1 = 28; % Length of LSCF, um
l2 = 72; % Length of GDC, um
Length=l1+l2; % um
int_wid=0;% Define the width of the interface region - 0 is the default value

delta_x=1; % Define the spatial step [um]
delta_t=1; % Define the time step [s]

steps_t=Duration*3600/delta_t+1; % Define number of time steps
steps_x=Length/delta_x+1; % Define number of spatial nodes

x=0:delta_x:Length; % Define the length of the profile in steps of delta_x
sigmaLSCF=DLSCF*delta_t/2*delta_x^2; % Calculate sigma of LSCF
sigmaGDC=DGDC*delta_t/2*delta_x^2; % Calculate sigma of GDC
sigma_int=D_int*delta_t/2*delta_x^2;% Calculate sigma of the interface
h=k/DLSCF; %calculate h

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set up a sub-matrix L which contains the information regarding the surface exchange as well as diffusion through the LSCF layer

L = full(gallery('tridiag',int8(l1/delta_x),sigmaLSCF,1-2*sigmaLSCF,sigmaLSCF)); % Bulk diffusion parameters
L(1,1) = 1-2*sigmaLSCF*(1+delta_x*h); %Surface boundary condition
L(1,2) = 2*sigmaLSCF; %Surface boundary condition


% Set up a sub-matrix E for the bulk diffusion properties of the
% electrolyte

E = full(gallery('tridiag',int8((l2-delta_x)/delta_x),sigmaGDC,1-2*sigmaGDC,sigmaGDC)); %Bulk diffusion parameters
E(end,end) = 1+sigmaGDC; %Mirror boundary condition parameters
E(end,end-1) = -2*sigmaGDC;
E(end,end-2) = sigmaGDC;

% Set up the main matrix which combines L, E and interface properties into one

An = full(gallery('tridiag',steps_x,sigmaGDC,1-2*sigmaGDC,sigmaGDC)); % Create a tridiagonal basis using the electrolyte values
An(1:l1/delta_x,1:l1/delta_x)=L; %integrate the sub-matrix L into the main matrix
An(Length-(l2-1)/delta_x+2:end,Length-(l2-1)/delta_x+2:end)=E; %integrate the sub-matrix E into the main matrix

%Set up the part of the matrix at the interface with varying width

for l=l1:l1+int_wid % Creates a loop which ensures that the tridiagonal matrix at the interface is set up with varying interface width
    An(l,l-1) = sigma_int;
    An(l,l) = 1-2*sigma_int;
    An(l,l+1) = sigma_int;
end

An(l1+int_wid,l1+1+int_wid) = sigmaGDC;
An(l1+1+int_wid,l1+int_wid) = sigma_int; % Manual modification to ensure the matrix is set up correctly


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Setting up the transformation matrix for the next timestep A(n+1), which
%is just the same as the previous matrix but with flipped signs.

Ln1 = full(gallery('tridiag',int8(l1/delta_x),-sigmaLSCF,1+2*sigmaLSCF,-sigmaLSCF));
Ln1(1,1) = 1+2*sigmaLSCF*(1+delta_x*h);
Ln1(1,2) = -2*sigmaLSCF;

En1 = full(gallery('tridiag',int8((l2-delta_x)/delta_x),-sigmaGDC,1+2*sigmaGDC,-sigmaGDC));
En1(end,end) = 1-sigmaGDC;
En1(end,end-1) = 2*sigmaGDC;
En1(end,end-2) = -sigmaGDC;

An1 = full(gallery('tridiag',steps_x,-sigmaGDC,1+2*sigmaGDC,-sigmaGDC));
An1(1:l1/delta_x,1:l1/delta_x)=Ln1;
An1(Length-(l2-1)/delta_x+2:end,Length-(l2-1)/delta_x+2:end)=En1;


for l=l1:l1+int_wid
    An1(l,l-1) = -sigma_int;
    An1(l,l) = 1+2*sigma_int;
    An1(l,l+1) = -sigma_int;
end


An1(l1+int_wid,l1+1+int_wid) = -sigmaGDC;
An1(l1+1+int_wid,l1+int_wid) = -sigma_int;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

G = zeros(steps_x, 1);
G(1) = 4*delta_x*sigmaLSCF*h; %surface boundary condition vector G

C = zeros(steps_x,1);%Initial position vector, 0 everywhere

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Run the actual iterative matrix calculation

for t=1:steps_t 
    
    C = An1\(An*C+G);
    plot(x,C,'x-');
    xlim([0 Length])
    ylim([0 inf]);
    xlabel('Distance, x / microns')
    ylabel('Concentration, C')
    drawnow
    
end


% Optional code to produce noisy data that can be used to test fitting
% functions

% snr = 50;
% out = awgn(C,snr);varargin;
%
% plot(out)
%
% writematrix(out,'GenNoiseInt.txt')
%