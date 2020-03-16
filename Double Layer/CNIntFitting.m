Dk0 = [9,264,100,0.017]; %Initial guesses for DLSCF, DGDC, DINT and k in that order


options = optimset('MaxFunEvals',2000, 'TolFun',1e-6,'TolX',1e-6);
[Dk, fminres] = fminsearch(@SUMSQInt,Dk0,options);

disp(['DLSCF = ',num2str(Dk(1)),'  DGDC = ',num2str(Dk(2)),' DINT = ',num2str(Dk(3)),' k = ',num2str(Dk(4))]);

cn=CNInt(Dk);
plot(cn)
hold on
GenNoise=load('DBLayerProfile.txt');
plot(GenNoise, '.')
legend('Fit','Data')
xlabel('Depth / um')
ylabel('Isotopic Fraction')