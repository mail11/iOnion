data=load('TripleLayerProfile.txt');

fun = @(Dk) sum((data-CNIntTL(Dk)).^2);

Dk0 = [0.01,0.01,0.1,0.01,1,0.001]; %Initial guesses for D1, D2, D3, DINT1, DINT 2 and k in that order


options = optimset('MaxIter',2000,'MaxFunEvals',3000, 'TolFun',1e-8,'TolX',1e-8);
[Dk, fminres] = fminsearch(fun,Dk0,options);

disp(['D1 = ',num2str(Dk(1)),'  D2 = ',num2str(Dk(2)),' D3 = ',num2str(Dk(3)),' DINT1 = ',num2str(Dk(4)),' DINT2 = ',num2str(Dk(5)),' k = ',num2str(Dk(6))]);

cn=CNIntTL(Dk);
plot(cn)
hold on
GenNoise=load('TripleLayerProfile.txt');
plot(GenNoise, '.')
legend('Fit','Data')
xlabel('Depth / um')
ylabel('Isotopic Fraction')