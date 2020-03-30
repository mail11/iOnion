function [Dk,data,cn] = CNIntFitting(C,data)

%% Define global variables
% This ensures that the following variables are the same across all
% function files

global steps_x;
global data;
global L_1;
global L_2;
global L_3;
global answer;
global Duration;
DeveloperMode = 1;

if DeveloperMode == 0
    %% Getting the data file
    
    prompt = {'Enter length of the profile in micrometers','Enter number of layers','Experiment time in h'}; % Creates a pop-up window that asks the user to enter the physical length of the diffusion profile in microns
    dlgtitle = 'Input';
    dims = 1;
    definput = {'1.0','2','0.083'};
    answer = inputdlg(prompt,dlgtitle,dims,definput);
    Duration = str2double(answer(3)); % Defines the experiment time
    
    [FileName, PathName] = uigetfile('*.*'); % Call a pop-up UI feature that allows the user to select data files
    data=transpose(textread([PathName,'\',FileName],'%n')'); % Reads the data file
    
    %% Extracting information from the data file
    
    if str2double(answer(2))== 1 % If there is only one single material
        
        scale = (length(data)-1)/str2double(answer(1));
        L_1 = length(data)/scale;
        L_2 = 0;
        L_3 = 0;
        
    elseif str2double(answer(2)) == 2 % Set the case for two layers and only picking one data point
        
        plot(data)
        [x] = ginput(1); % Asks the user to select a point on the graph to mark the interface location, then extracts the x value of that point
        scale = (length(data)-1)/str2double(answer(1)); % Sets a scaling factor based on the user input to allow conversion of the x values from the datafile into real lengths
        L_1 = x(1)/scale; % Determines the length of the first layer based on user input
        L_2 = (length(data)-1)/scale-L_1; % Determines the length of the second layer
        L_3 = 0;
        close
        
    elseif str2double(answer(2)) == 3 % Set the case for three layers and picking two data points
        
        plot(data)
        [x] = ginput(2); % Asks the user to select a point on the graph to mark the interface location, then extracts the x value of that point
        scale = (length(data)-1)/str2double(answer(1));
        L_1 = x(1)/scale;
        L_2 = x(2)/scale;
        L_3 = (length(data)-1)/scale-(L_1+L_2);
        close
    end
    
    %% Setting up and running fminsearch
    
    fun = @(Dk) sum((data-CNInt(Dk)).^2); % set up the sum squared function
    
    if str2double(answer(2)) == 1 % Condition for a single material
        Dk0 = log([0.1,0.1]); % D and k
    elseif    str2double(answer(2)) == 2 % Ensures that 4 variables are run in fminsearch
        Dk0 = log([0.009,0.264,0.00100,0.0017]); % Initial guesses for DLSCF, DGDC, DINT and k in that order
    elseif str2double(answer(2)) == 3 % Ensures that 6 variables are run
        Dk0 = log([0.1,0.1,0.1,0.1,1,0.1]); %Initial guesses for D1, D2, D3, DINT1, DINT 2 and k in that order
    end
    
    options = optimset('MaxIter',2000,'MaxFunEvals',3000, 'TolFun',1e-8,'TolX',1e-8); % set options for fminsearch
    [Dk, fminres] = fminsearch(fun,Dk0,options); % the actual fminsearch
    
    % Change the output depending on how many variables are fitted
    
    if str2double(answer(2)) == 1
        disp(['D1 = ',num2str((10^(Dk(1)))),' k = ',num2str((10^(Dk(2))))]);
    elseif str2double(answer(2)) == 2
        disp(['D1 = ',num2str((10^(Dk(1)))),'  D2 = ',num2str((10^(Dk(2)))),' DINT = ',num2str((10^(Dk(3)))),' k = ',num2str((10^(Dk(4))))]); % display the output variables
    elseif str2double(answer(2)) == 3
        disp(['D1 = ',num2str((10^(Dk(1)))),'  D2 = ',num2str((10^(Dk(2)))),' D3 = ',num2str((10^(Dk(3)))),' DINT1 = ',num2str((10^(Dk(4)))),' DINT2 = ',num2str((10^(Dk(5)))),' k = ',num2str((10^(Dk(6))))]); % display the output variables
    end
    
    writematrix(Dk,'FittedValues.xlsx') % Write the values out onto an Excel file which can then be further processed
    
elseif DeveloperMode == 1
    
    data=load('GenNoiseInt.txt');
    L_1 = 0.3;
    L_2 = 0.7;
    L_3 = 0;
    answer = cell(3,1);
    answer(2) = {'2'};
    str2double(answer(2));
    Duration = 0.083;
    
    fun = @(Dk) sum((data-CNInt(Dk)).^2);
    
    Dk0 = log([0.009,0.264,0.00100,0.0017]); %Initial guesses for DLSCF, DGDC, DINT and k in that order
    
    
    options = optimset('MaxFunEvals',2000, 'TolFun',1e-6,'TolX',1e-6);
    [Dk, fminres] = fminsearch(fun,Dk0,options);
    

    disp(['DLSCF = ',num2str((10^(Dk(1)))),'  DGDC = ',num2str((10^(Dk(2)))),' DINT = ',num2str((10^(Dk(3)))),' k = ',num2str((10^(Dk(4)))));
    
end
%% plotting

cn=CNInt(Dk); % Calls the final output from the fitting algorithm
x=linspace(0,L_1+L_2+L_3,steps_x); % Set up x values in accordance with the real scale of the data
plot(x,cn)
hold on
plot(x,data, '.')
legend('Fit','Data')
xlabel('Depth / um')
ylabel('Isotopic Fraction')

end