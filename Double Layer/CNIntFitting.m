%% Getting the data file

[FileName, PathName] = uigetfile('*.*');
global data;
data=transpose(textread([PathName,'\',FileName],'%n')');

%% Setting up fminsearch

fun = @(Dk) sum((data-CNInt(Dk)).^2); % set up the sum squared function

Dk0 = [0.009,0.264,0.00100,0.0017]; % Initial guesses for DLSCF, DGDC, DINT and k in that order


options = optimset('MaxIter',2000,'MaxFunEvals',3000, 'TolFun',1e-8,'TolX',1e-8); % set options for fminsearch
[Dk, fminres] = fminsearch(fun,Dk0,options); % the actual fminsearch

disp(['DLSCF = ',num2str(Dk(1)),'  DGDC = ',num2str(Dk(2)),' DINT = ',num2str(Dk(3)),' k = ',num2str(Dk(4))]); % display the output variables

writematrix(Dk,'FittedValues.xlsx') % Write the values out onto an Excel file which can then be further processed

%% plotting

cn=CNInt(Dk);
plot(cn)
hold on
plot(data, '.')
legend('Fit','Data')
xlabel('Depth / um')
ylabel('Isotopic Fraction')