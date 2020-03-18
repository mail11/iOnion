%% Getting the data file

[FileName, PathName] = uigetfile('*.*');
global data;
data=transpose(textread([PathName,'\',FileName],'%n')');

%% Setting up fminsearch

fun = @(Dk) sum((data-CNIntTL(Dk)).^2); % set up the sum squared function

Dk0 = [0.01,0.01,0.1,0.01,1,0.001]; %Initial guesses for D1, D2, D3, DINT1, DINT 2 and k in that order


options = optimset('MaxIter',2000,'MaxFunEvals',3000, 'TolFun',1e-8,'TolX',1e-8); % set options for fminsearch
[Dk, fminres] = fminsearch(fun,Dk0,options); % the actual fminsearch

disp(['D1 = ',num2str(Dk(1)),'  D2 = ',num2str(Dk(2)),' D3 = ',num2str(Dk(3)),' DINT1 = ',num2str(Dk(4)),' DINT2 = ',num2str(Dk(5)),' k = ',num2str(Dk(6))]); % display the output variables


writematrix(Dk,'FittedValues.xlsx') % Write the values out onto an Excel file which can then be further processed

%% Plotting

cn=CNIntTL(Dk);
plot(cn)
hold on
GenNoise=load('TripleLayerProfile.txt');
plot(GenNoise, '.')
legend('Fit','Data')
xlabel('Depth / um')
ylabel('Isotopic Fraction')