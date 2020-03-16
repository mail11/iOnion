function C = CNInt(Dk)

% DLSCF=10; % Define diffusivity of material [um^2/s] (LSCF)
% DGDC=50;
% D_int=10;
% k=0.0469; %define k [um/s]
Duration=.083; %hours
l1 = 25; %length of LSCF, um
l2 = 75; %length of GDC, um
Length=l1+l2; %um
int_wid=2;

delta_x=2; % Define the spatial step [um]
delta_t=2; % Define the time step [s] 

steps_t=Duration*3600/delta_t+1; % Define number of time steps
steps_x=Length/delta_x+1; % Define number of spatial nodes

x=0:delta_x:Length;
sigmaLSCF=Dk(1)*delta_t/2*delta_x^2; %Calculate sigma
sigmaGDC=Dk(2)*delta_t/2*delta_x^2; %Calculate sigma
sigma_int=Dk(3)*delta_t/2*delta_x^2;%Interfacial Resistance
h=Dk(4)/Dk(1); %calculate h

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
L = full(gallery('tridiag',l1/delta_x,sigmaLSCF,1-2*sigmaLSCF,sigmaLSCF));
L(1,1) = 1-2*sigmaLSCF*(1+delta_x*h);
L(1,2) = 2*sigmaLSCF;


E = full(gallery('tridiag',(l2-1)/delta_x,sigmaGDC,1-2*sigmaGDC,sigmaGDC));
E(end,end) = 1+sigmaGDC;
E(end,end-1) = -2*sigmaGDC;
E(end,end-2) = sigmaGDC;

An = full(gallery('tridiag',steps_x,sigmaGDC,1-2*sigmaGDC,sigmaGDC));
An(1:l1/delta_x,1:l1/delta_x)=L;



An(Length-(l2-1)/delta_x+2:end,Length-(l2-1)/delta_x+2:end)=E;



for l=l1:l1+int_wid
    An(l,l-1) = sigma_int;
    An(l,l) = 1-2*sigma_int;
    An(l,l+1) = sigma_int;
end


% An(l1+1,l1) = sigmaLSCF;
 An(l1+int_wid,l1+1+int_wid) = sigmaGDC;
 An(l1+1+int_wid,l1+int_wid) = sigma_int;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Ln1 = full(gallery('tridiag',l1/delta_x,-sigmaLSCF,1+2*sigmaLSCF,-sigmaLSCF));
Ln1(1,1) = 1+2*sigmaLSCF*(1+delta_x*h);
Ln1(1,2) = -2*sigmaLSCF;


En1 = full(gallery('tridiag',(l2-1)/delta_x,-sigmaGDC,1+2*sigmaGDC,-sigmaGDC));
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

% An1(l1+1,l1) = -sigmaLSCF;
An1(l1+int_wid,l1+1+int_wid) = -sigmaGDC;
An1(l1+1+int_wid,l1+int_wid) = -sigma_int;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

G = zeros(steps_x, 1);
G(1) = 4*delta_x*sigmaLSCF*h; %surface boundary condition vector G

C = zeros(steps_x,1);%Initial position vector, 0 everywhere

for t=1:steps_t
    
    C = An1\(An*C+G);
    plot(x,C);
    xlim([0 Length])
    ylim([0.06 inf]);
    xlabel('Distance, x / microns')
    ylabel('Concentration, C')
%     drawnow
    
 
end
    

end