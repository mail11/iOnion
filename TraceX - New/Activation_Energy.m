D_LSCF = transp([9.58E-18;9.96E-17;7.64E-16;4.58E-15;2.24E-14;9.20E-14;3.27E-13;
    1.03E-12; 2.90E-12]);
D_GDC = transp([1.59E-14; 4.59E-14; 1.14E-13; 2.61E-13; 5E-13; 1E-12; 2E-12; 3E-12;
    5E-12]);
k_LSCF = transp([8.40E-11; 1.76E-10; 3.54E-10;  6.05E-10; 1.17E-09; 4.69E-09; 1.42E-08
    4.45E-08; 1.14E-07]);
kB = 8.6173e-05;
T = (600:50:800)+273.5;
Tinv = 1./T;
y = polyfit(Tinv,log(k_LSCF),1); 
D0 = exp(y(2));
EA = -y(1)*kB;
disp(['D0 = ',num2str(D0),' EA (eV) = ',num2str(EA)]);

%D = D0*exp(-EA/k*T));
% D0*exp((-EA/kB)*Tinv))
% 0.0049255*exp((-1.9659/kB)*Tinv))



