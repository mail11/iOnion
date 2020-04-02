%This script produces a Crank-Nicolson simulation of a stack of two layers
%of material. The main purpose is to simulate a diffusion profile of a system with a GDC
%electrolyte and an LSCF cathode as well as an interface between the two.

clear;close all

%% Define parameters
% Material properties
D_1=.4; % Define diffusivity of cathode material [um^2/s] (LSCF)
D_2=4;  % Define diffusivity of electrolyte material [um^2/s] (GDC)
k=0.0772; % Define surface exchange coefficient of LSCF [um/s]
D_int=.1; %Define interfacial resistance [s/um]
int_width=0;% Define the width of the interface region (0 is the default value)

% Experimental setup
Duration=.2; % Time of exchange in hours (found with Kiloran correction)
L = [20;80]; % Vector of the layer lengths, first element is the first layer

% Simulation parameters
delta_x=1; % Define the spatial step [um]
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

sub = zeros(round(steps_x)-1,1);
sub(1:L(1)/delta_x-1,1)=sigma_1;
sub(L(1)/delta_x:L(1)/delta_x+1,1)=sigma_int;
sub(L(1)/delta_x+1:end,1)=sigma_2;

sup = zeros(round(steps_x)-1,1);
sup(1:L(1)/delta_x-2,1)=sigma_1;
sup(L(1)/delta_x-1:L(1)/delta_x,1)=sigma_int;
sup(L(1)/delta_x+0:end,1)=sigma_2;

s1 = zeros(round(steps_x),1);
s1(1:L(1)/delta_x-1,1)=1-2*sigma_1;
s1(L(1)/delta_x:L(1)/delta_x+1,1)=1-2*sigma_int;
s1(L(1)/delta_x+1:end,1)=1-2*sigma_2;

A = full(gallery('tridiag',sub,s1,sup));
A(1,1) = 1-2*sigma_1*(1+delta_x*h); %Surface boundary condition
A(1,2) = 2*sigma_1; %Surface boundary condition
A(end,end) = 1+sigma_2; %Mirror boundary condition parameters
A(end,end-1) = -2*sigma_2;
A(end,end-2) = sigma_2;

A_n = full(gallery('tridiag',-sub,-s1+2,-sup));
A_n(1,1) = 1+2*sigma_1*(1+delta_x*h); %Surface boundary condition
A_n(1,2) = -2*sigma_1; %Surface boundary condition
A_n(end,end) = 1-sigma_2; %Mirror boundary condition parameters
A_n(end,end-1) = 2*sigma_2;
A_n(end,end-2) = -sigma_2;

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
C_Di = zeros(round(steps_x),1);


%% Run the iterative matrix calculation

for t=1:steps_t
    C = A_n\(A*C+G);
    C_Di = A_n_Di\(A_Di*C_Di+G);
    plot(x,C_Di,x,C,x,(C+C_Di)/2);
    xlim([0 Length])
    ylim([-0 inf]);
    xlabel('Distance, x / microns')
    ylabel('Concentration, C')
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

%
% Layer1 = full(gallery('tridiag',round(L(1)/delta_x),sigma_1,1-2*sigma_1,sigma_1)); % Bulk diffusion parameters
% Layer1(1,1) = 1-2*sigma_1*(1+delta_x*h); %Surface boundary condition
% Layer1(1,2) = 2*sigma_1; %Surface boundary condition
%
% % Set up a sub-matrix E for the bulk diffusion properties of the electrolyte
% Layer2 = full(gallery('tridiag',round((L(2)-delta_x)/delta_x)+2,sigma_2,1-2*sigma_2,sigma_2)); %Bulk diffusion parameters
% Layer2(end,end) = 1+sigma_2; %Mirror boundary condition parameters
% Layer2(end,end-1) = -2*sigma_2;
% Layer2(end,end-2) = sigma_2;
%
% % Set up the main matrix which combines L, E and interface properties into one
% A = full(gallery('tridiag',round(steps_x),sigma_2,1-2*sigma_2,sigma_2)); % Create a tridiagonal basis using the electrolyte values
% A(1:L(1)/delta_x,1:L(1)/delta_x)=Layer1; %integrate the sub-matrix L into the main matrix
% A(round((Length-L(2))/delta_x)+1:end,round((Length-L(2))/delta_x)+1:end)=Layer2; %integrate the sub-matrix E into the main matrix
%
% % %Set up the part of the matrix at the interface with varying width
% for l=L(1)/delta_x:L(1)/delta_x+int_width % Creates a loop which ensures that the tridiagonal matrix at the interface is set up with varying interface width
%     A(l,l-1) = sigma_int;
%     A(l,l) = 1-2*sigma_int;
%     A(l,l+1) = sigma_int;
% end
%
%
% A(L(1)/delta_x+int_width,L(1)/delta_x+1+int_width) = sigma_2;
% A(L(1)/delta_x+1+int_width,L(1)/delta_x+int_width) = sigma_int; % Manual modification to ensure the matrix is set up correctly


% A(L(1)/delta_x,L(1)/delta_x-1)=sigma_1+0.25*(sigma_int-sigma_1);
% A(L(1)/delta_x,L(1)/delta_x)=1-2*sigma_1;
% A(L(1)/delta_x,L(1)/delta_x+1)=sigma_1-0.25*(sigma_int-sigma_1);
%
% A(L(1)/delta_x+1,L(1)/delta_x)=sigma_int+0.25*(sigma_2-sigma_1);
% A(L(1)/delta_x+1,L(1)/delta_x+1)=1-2*sigma_int;
% A(L(1)/delta_x+1,L(1)/delta_x+2)=sigma_int-0.25*(sigma_2-sigma_1);
%
% A(L(1)/delta_x+2,L(1)/delta_x+1)=sigma_2+0.25*(sigma_2-sigma_int);
% A(L(1)/delta_x+2,L(1)/delta_x+2)=1-2*sigma_2;
% A(L(1)/delta_x+2,L(1)/delta_x+3)=sigma_2-0.25*(sigma_2-sigma_int);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Setting up the transformation matrix for the next timestep A(n+1), which
%is just the same as the previous matrix but with flipped signs.

% Layer1_n = full(gallery('tridiag',round(L(1)/delta_x),-sigma_1,1+2*sigma_1,-sigma_1));
% Layer1_n(1,1) = 1+2*sigma_1*(1+delta_x*h);
% Layer1_n(1,2) = -2*sigma_1;
%
% Layer2_n = full(gallery('tridiag',round((L(2)-delta_x)/delta_x)+2,-sigma_2,1+2*sigma_2,-sigma_2));
% Layer2_n(end,end) = 1-sigma_2;
% Layer2_n(end,end-1) = 2*sigma_2;
% Layer2_n(end,end-2) = -sigma_2;
%
% A_n = full(gallery('tridiag',round(steps_x),-sigma_2,1+2*sigma_2,-sigma_2));
% A_n(1:L(1)/delta_x,1:L(1)/delta_x)=Layer1_n;
% A_n(round((Length-L(2))/delta_x)+1:end,round((Length-L(2))/delta_x)+1:end)=Layer2_n;
%
% for l=L(1)/delta_x:(L(1)+int_width)/delta_x
%     A_n(l,l-1) = -sigma_int;
%     A_n(l,l) = 1+2*sigma_int;
%     A_n(l,l+1) = -sigma_int;
% end
%
% A_n(L(1)/delta_x+int_width,L(1)/delta_x+1+int_width) = -sigma_2;
% A_n(L(1)/delta_x+1+int_width,L(1)/delta_x+int_width) = -sigma_int;



% A_n(L(1)/delta_x,L(1)/delta_x-1)=-sigma_1+0.25*(sigma_int-sigma_1);
% A_n(L(1)/delta_x,L(1)/delta_x)=1+2*sigma_1;
% A_n(L(1)/delta_x,L(1)/delta_x+1)=-sigma_1-0.25*(sigma_int-sigma_1);
%
% A_n(L(1)/delta_x+1,L(1)/delta_x)=-sigma_int+0.25*(sigma_2-sigma_1);
% A_n(L(1)/delta_x+1,L(1)/delta_x+1)=1+2*sigma_int;
% A_n(L(1)/delta_x+1,L(1)/delta_x+2)=-sigma_int-0.25*(sigma_2-sigma_1);
%
% A_n(L(1)/delta_x+2,L(1)/delta_x+1)=-sigma_2+0.25*(sigma_2-sigma_int);
% A_n(L(1)/delta_x+2,L(1)/delta_x+2)=1+2*sigma_2;
% A_n(L(1)/delta_x+2,L(1)/delta_x+3)=-sigma_2-0.25*(sigma_2-sigma_int);