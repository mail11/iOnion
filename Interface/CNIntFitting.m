function [r0,Dk,data,cn] = CNIntFitting(Duration,data,Length_um,L) 
%% Define global variables
% This ensures that the following variables are the same across all
% function files

global steps_x;
global delta_x;
global data;
global Duration;
global Layers;
global Length_um;
global L;
L = [0;0;0];
DeveloperMode = 1;

if DeveloperMode == 0
    %% Getting the data file
    
    prompt = {'Enter length of the profile in micrometers','Enter number of layers','Experiment time in h','Gas concentration','Background concentration'}; % Creates a pop-up window that asks the user to enter the physical length of the diffusion profile in microns
    dlgtitle = 'Input';
    dims = 1;
    definput = {'1.0','2','0.083','0.98','0.02'};
    answer = inputdlg(prompt,dlgtitle,dims,definput);
    Length_um=str2double(answer(1)); % Provides the physical scaling of the fitting data
    Layers=str2double(answer(2)); % Determines for how many layers the optimisation should run
    Duration=str2double(answer(3));% Defines the experiment time
    C_gas=str2double(answer(4)); % Defines the concentration of O18 in the gas
    C_bg=str2double(answer(5)); % Defines the background concentration of O18 in the native oxide lattice
    
    [FileName, PathName] = uigetfile('*.*'); % Call a pop-up UI feature that allows the user to select data files
    data=(transpose(textread([PathName,'\',FileName],'%n')')-C_bg)/(C_gas-C_bg); % Reads the data file
    
    %% Extracting information from the data file
    
    if Layers == 1 % If there is only one single material
        
        scale = (length(data)-1)/Length_um;
        L(1) = length(data)/scale;
        L(2) = 0;
        L(3) = 0;
        
    elseif Layers == 2 % Set the case for two layers and only picking one data point
        
        plot(data)
        [x] = ginput(1); % Asks the user to select a point on the graph to mark the interface location, then extracts the x value of that point
        scale = (length(data)-1)/Length_um; % Sets a scaling factor based on the user input to allow conversion of the x values from the datafile into real lengths
        L(1) = x(1)/scale; % Determines the length of the first layer based on user input
        L(2) = (length(data)-1)/scale-L(1); % Determines the length of the second layer
        L(3) = 0;
        close
        
    elseif Layers == 3 % Set the case for three layers and picking two data points
        
        plot(data)
        [x] = ginput(2); % Asks the user to select a point on the graph to mark the interface location, then extracts the x value of that point
        scale = (length(data)-1)/Length_um;
        L(1) = x(1)/scale;
        L(2) = x(2)/scale;
        L(3) = (length(data)-1)/scale-(L(1)+L(2));
        close
    end
    
    %% Setting up and running fminsearch
    
    fun = @(Dk) sum((data-CNInt(Dk)).^2); % set up the sum squared function
    
    if Layers == 1 % Condition for a single material
        Dk0 = log([0.1,0.1]); % D and k
    elseif    Layers == 2 % Ensures that 4 variables are run in fminsearch
        Dk0 = log([0.009,0.264,0.00100,0.0017]); % Initial guesses for DLSCF, DGDC, DINT and k in that order
    elseif Layers == 3 % Ensures that 6 variables are run
        Dk0 = log([0.1,0.1,0.1,0.1,0.1,0.01]); %Initial guesses for D1, D2, D3, DINT1, DINT 2 and k in that order
    end
    
    options = optimset('MaxIter',2000,'MaxFunEvals',3000, 'TolFun',1e-8,'TolX',1e-8); % set options for fminsearch
    [Dk, fminres] = fminsearch(fun,Dk0,options); % the actual fminsearch
    
    % Change the output depending on how many variables are fitted
    
    if Layers == 1
        disp(['D1 = ',num2str((10^(Dk(1)))),' k = ',num2str((10^(Dk(2))))]);
    elseif Layers == 2
        disp(['D1 = ',num2str((10^(Dk(1)))),'  D2 = ',num2str((10^(Dk(2)))),' DINT = ',num2str((10^(Dk(3)))),' k = ',num2str((10^(Dk(4))))]); % display the output variables
        r = (2/10^Dk(3)-1/10^Dk(1)-1/10^Dk(2))*delta_x
    elseif Layers == 3
        disp(['D1 = ',num2str((10^(Dk(1)))),'  D2 = ',num2str((10^(Dk(2)))),' D3 = ',num2str((10^(Dk(3)))),' DINT1 = ',num2str((10^(Dk(4)))),' DINT2 = ',num2str((10^(Dk(5)))),' k = ',num2str((10^(Dk(6))))]); % display the output variables
        r1 = (2/10^Dk(4)-1/10^Dk(1)-1/10^Dk(2))*delta_x
        r2 = (2/10^Dk(5)-1/10^Dk(2)-1/10^Dk(3))*delta_x
    end
    
    writematrix(Dk,'FittedValues.xlsx') % Write the values out onto an Excel file which can then be further processed
    
elseif DeveloperMode == 1
    
    
    C_bg = 0.02;
    C_gas = 0.98;
    data=(load('GenNoiseInt.txt')-C_bg)/(C_gas-C_bg);
    L(1) = 0.3;
    L(2) = 0.7;
    L(3) = 0;
    Layers = 2;
    Duration = 0.083;
    
    fun = @(Dk) sum((data-CNInt(Dk)).^2);
    
    Dk0 = log([0.1,0.1,0.1,0.1]); %Initial guesses for DLSCF, DGDC, DINT and k in that order
    
    
    options = optimset('MaxFunEvals',2000, 'TolFun',1e-6,'TolX',1e-6);
    [Dk, fminres] = fminsearch(fun,Dk0,options);
    
    
    disp(['DLSCF = ',num2str((10^(Dk(1)))),'  DGDC = ',num2str((10^(Dk(2)))),' DINT = ',num2str((10^(Dk(3)))),' k = ',num2str((10^(Dk(4))))]);
    r = (2/10^Dk(3)-1/10^Dk(1)-1/10^Dk(2))*delta_x
end
%% plotting

cn=CNInt(Dk); % Calls the final output from the fitting algorithm
x=linspace(0,L(1)+L(2)+L(3),steps_x); % Set up x values in accordance with the real scale of the data
plot(x,cn)
hold on
plot(x,data, '.')
legend('Fit','Data')
xlabel('Depth / um')
ylabel('Isotopic Fraction')

end