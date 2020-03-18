data=load('DBLayerProfile.txt');

fun = @(Dk) sum((data-CNInt(Dk)).^2);

Dk0 = [0.009,0.264,0.00100,0.0017]; %Initial guesses for DLSCF, DGDC, DINT and k in that order


options = optimset('MaxFunEvals',2000, 'TolFun',1e-6,'TolX',1e-6);
[Dk, fminres] = fminsearch(fun,Dk0,options);

disp(['DLSCF = ',num2str(Dk(1)),'  DGDC = ',num2str(Dk(2)),' DINT = ',num2str(Dk(3)),' k = ',num2str(Dk(4))]);

%% plotting

cn=CNInt(Dk);
plot(cn)
hold on
plot(data, '.')
legend('Fit','Data')
xlabel('Depth / um')
ylabel('Isotopic Fraction')