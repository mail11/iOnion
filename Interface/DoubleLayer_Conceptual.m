%This script produces a Crank-Nicolson simulation of a stack of two layers
%of material. The main purpose is to simulate a diffusion profile of a system with a GDC
%electrolyte and an LSCF cathode as well as an interface between the two.

clear;close all



%% Define parameters
% Material properties
D_1=4.58*1e-3; % Define diffusivity of cathode material [um^2/s] (LSCF)
D_2=2.61*1e-1;  % Define diffusivity of electrolyte material [um^2/s] (GDC)
k=6.05*1e-2; % Define surface exchange coefficient of LSCF [um/s]

for i = 1:11; r=2^i;

    %r = 10000000; %Define interfacial resistance [s/um]


int_width=0;% Define the width of the interface region (0 is the default value)

% Experimental setup
Duration=.00965; % Time of exchange in hours (found with Kiloran correction)
L = [3;10]; % Vector of the layer lengths, first element is the first layer

% Simulation parameters
delta_x=0.1; % Define the spatial step [um]
delta_t=2; % Define the time step [s] (this should be a derived parameter based on max value of sigma)

D_int=2/(2*r/delta_x+1/D_1+1/D_2); 

% Derived properties
%D_int=2*(r/delta_x+1/D_1+1/D_2)^-1; % Define the diffusivity of the interface [um^2/s]
Length=L(1)+L(2); % Total sample length [um]
steps_x=Length/delta_x+1; % Define number of spatial nodes
steps_t=Duration*3600/delta_t+1; % Define number of time steps
x=0:delta_x:Length; % Define the length of the profile in steps of delta_x
x_n = x/(2*sqrt(D_1*Duration*3600));
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

A_nI = inv(A_n);
A_nI_Di = inv(A_n_Di);

for t=1:steps_t
    C = A_nI*(A*C+G);
    C_Di = A_nI_Di*(A_Di*C_Di+G);
end


%plot(x,C_Di,x,C,x,(C+C_Di)/2); % Plot both boundary conditions as well as their mean
plot(x_n,(C+C_Di)/2,'Marker','none','LineWidth',1,'Color','[0.65 0.66 0.89]') % Plot only the mean
xlim([0 10])
ylim([-0 inf]);
set(gca,'XTick',(0:0.2:2.5),'FontSize',24)
set(gcf,'color','w');
xlabel('Normalised Profile Depth, {\it x''}')
ylabel('Isotopic Fraction, {\it C}')
% dim = [.53 .4 .25 .2];
% str = {'k = 2.92*10^{-10} m/s ','r = 3.63*10^{-4} s/m','Exchange Time = 0.2 hours'};
% annotation('textbox',dim,'String',str,'FitBoxToText','on')


% yyaxis right
% %semilogy(0.5*delta_x+x(1,1:L(1)/delta_x),10^-12*(sub(1:L(1)/delta_x)*2*delta_x^2)/delta_t,'LineStyle','-','LineWidth',2,'Color','[0.733 0.52 0.67]','Marker','none','MarkerFaceColor','[0.733 0.52 0.67]','MarkerEdgeColor','[0.733 0.52 0.67]','MarkerSize',4);
% %hold on
% semilogy(0.5*delta_x+x(1,L(1)/delta_x+1:100),10^-12*(sub(L(1)/delta_x+1:100)*2*delta_x^2)/delta_t,'LineStyle','-','LineWidth',1.5,'Color','[0.733 0.52 0.67]','Marker','none','MarkerFaceColor','[0.733 0.52 0.67]','MarkerEdgeColor','[0.733 0.52 0.67]','MarkerSize',4);

% ylabel ('Tracer Diffusion Coefficient {\it D*}/m^2/s','Color','k')
% set(gca,'ycolor','k')
%set(gca,'YTick',[])
%set(gca,'YTickLabel')

% xline((L(1)-delta_x)/(2*D_1*3600*Duration),'--','Color','k');
% xline(L(1),'--','Color','k');
y = 10^-12*D_1;
% line([0,L(1)-delta_x],[y,y],'LineWidth',1.5,'Color','[0.733 0.52 0.67]')
hold on
% y = 10^-12*D_int;
% line([L(1)-delta_x,L(1)],[y,y],'LineWidth',1.5,'Color','[0.733 0.52 0.67]')
% hold on
% y = 10^-12*D_2;
% line([L(1),L(1)+L(2)],[y,y],'LineWidth',1.5,'Color','[0.733 0.52 0.67]')


% legend('Simulated Profile','Interface')

end

%r = (2/D_int-1/D_1-1/D_2)*delta_x/2

%hold on
% drawnow
% 

% 
% %Optional code to produce noisy data that can be used to test fitting
% %functions
% 
% % snr = 50;
% % out = awgn(C,snr);varargin;
% %
% % plot(out)
% %
% % writematrix(out,'GenNoiseInt.txt')



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
