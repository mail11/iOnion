%This script produces a Crank-Nicolson simulation of a stack of two layers
%of material. The main purpose is to simulate a diffusion profile of a system with a GDC
%electrolyte and an LSCF cathode as well as an interface between the two.

clear;close all



%% Define parameters
% Material properties
D_1=0.0000138; % Define diffusivity of cathode material [um^2/s] (LSCF)
D_2=0.0000138;  % Define diffusivity of electrolyte material [um^2/s] (GDC)
k = 0.0001;
for i = 1:25; r=2^i;
%r = 0; %Define interfacial resistance [s/um]
int_width=0;% Define the width of the interface region (0 is the default value)

% Experimental setup

%for i = 0:15; Duration=i*0.5;
Duration=1.5; % Time of exchange in hours (found with Kiloran correction)
L = [0.5;0.5]; % Vector of the layer lengths, first element is the first layer

% Simulation parameters
delta_x=0.005; % Define the spatial step [um]
delta_t=2; % Define the time step [s] (this should be a derived parameter based on max value of sigma)

D_int=2/(2*r/delta_x+1/D_1+1/D_2); 

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
A(1,1) = sigma_2;
A(1,2) = -2*sigma_2;
A(1,3) = 1+sigma_2; % Left mirror boundary condition parameters
A(end,end) = 1+sigma_2; % Right mirror boundary condition parameters
A(end,end-1) = -2*sigma_2;
A(end,end-2) = sigma_2;

% Create the matrix at the future timestep
A_n = full(gallery('tridiag',-sub,s1_n,-sup));
A_n(1,1) = -sigma_2;
A_n(1,2) = 2*sigma_2;
A_n(1,3) = 1-sigma_2; 
A_n(end,end) = 1-sigma_2; %Mirror boundary condition parameters
A_n(end,end-1) = 2*sigma_2;
A_n(end,end-2) = -sigma_2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

C = zeros(round(steps_x),1);%Initial position vector, 0 everywhere
C(1:L(1)/delta_x) = 1;
%C_Di = zeros(round(steps_x),1);% Initial position vector for the Dirichlet bc


%% Run the iterative matrix calculation


for t=1:steps_t
    C = A_n\(A*C);
    %C_Di = A_n_Di\(A_Di*C_Di+G);
    %drawnow
end


%plot(x,C_Di,x,C,x,(C+C_Di)/2); % Plot both boundary conditions as well as their mean
plot(x,C,'Marker','none','LineWidth',1,'Color','[0.65 0.66 0.89]') % Plot only the mean
xlim([0 Length])
ylim([-0 inf]);
set(gca,'XTick',(0:Length/10:Length))
xlabel('Profile Depth, {\it x} / microns')
ylabel('Isotopic Fraction, {\it C}')
% dim = [.53 .4 .25 .2];
% str = {'k = 2.92*10^{-10} m/s ','r = 3.63*10^{-4} s/m','Exchange Time = 0.2 hours'};
% annotation('textbox',dim,'String',str,'FitBoxToText','on')


% yyaxis right
% %semilogy(0.5*delta_x+x(1,1:L(1)/delta_x),10^-12*(sub(1:L(1)/delta_x)*2*delta_x^2)/delta_t,'LineStyle','-','LineWidth',2,'Color','[0.733 0.52 0.67]','Marker','none','MarkerFaceColor','[0.733 0.52 0.67]','MarkerEdgeColor','[0.733 0.52 0.67]','MarkerSize',4);
% %hold on
% semilogy(0.5*delta_x+x(1,L(1)/delta_x+1:100),10^-12*(sub(L(1)/delta_x+1:100)*2*delta_x^2)/delta_t,'LineStyle','-','LineWidth',1.5,'Color','[0.733 0.52 0.67]','Marker','none','MarkerFaceColor','[0.733 0.52 0.67]','MarkerEdgeColor','[0.733 0.52 0.67]','MarkerSize',4);

%ylabel ('Tracer Diffusion Coefficient {\it D*}/m^2/s','Color','k')
set(gca,'ycolor','k')
%set(gca,'YTick',[])
%set(gca,'YTickLabel')

xline(L(1)-delta_x,'--','Color','k');
xline(L(1),'--','Color','k');
y = 10^-12*D_1;
% line([0,L(1)-delta_x],[y,y],'LineWidth',1.5,'Color','[0.733 0.52 0.67]')
hold on
% y = 10^-12*D_int;
% line([L(1)-delta_x,L(1)],[y,y],'LineWidth',1.5,'Color','[0.733 0.52 0.67]')
% hold on
% y = 10^-12*D_2;
% line([L(1),L(1)+L(2)],[y,y],'LineWidth',1.5,'Color','[0.733 0.52 0.67]')


%legend('Simulated Profile','D*','Interface')

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