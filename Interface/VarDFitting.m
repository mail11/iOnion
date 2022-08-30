function [Dk,data,cn] = VarDFitting(Duration,data,L)
%% Define global variables
global data;
global L;
global delta_x;
global delta_t;
global steps_x;
global Duration;
global Layers;
global C_gas;
global CM;
global xq;

%% Getting the data file

load CM
load xM

data = transpose(CM);
%% Extracting information from the data file


%% Setting up and running fminsearch

fun = @(Dk) sum((data-VarD(Dk)).^2); % set up the sum squared function

Dk0 = 1E-3*[1.0,0.00001,58,-0.0,-0.0]; %Initial guesses for D1, D3, k, A, B in that order
options = optimset('Display','iter','PlotFcns',@optimplotfval,'MaxFunEvals',2000,'TolFun',1e-4,'TolX',1e-4);

[Dk, fminres,exitflag,output] = fminsearch(fun,Dk0,options);

%% plotting

cn=VarD(Dk); % Calls the final output from the fitting algorithm
plot(xq,cn)
hold on
plot(xq,data, '.')
legend('Fit','Data')
xlabel('Depth / um')
ylabel('Isotopic Fraction')

disp(['D1 = ',num2str(Dk(1)),'  D2 = ',num2str(Dk(2)),' k = ',num2str(Dk(3)),' A = ',num2str(((1+exp(-Dk(4)))^-1)),' B = ',num2str((exp(Dk(5))))]);

end


