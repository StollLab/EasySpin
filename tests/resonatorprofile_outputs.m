function ok = test()

% Compare resonatorprofile() outputs with mathematical expressions
nu = linspace(9.2,9.8,101); % GHz
nu0 = 9.5; % GHz
Qu = 2000;
beta = 1;

% Transfer function
H0 = 1./(1+1i*Qu*(nu/nu0-nu0./nu));
H = resonatorprofile(nu,nu0,Qu,'transferfunction');

ok(1) = areequal(H0,H,1e-12,'abs');

% Voltage reflection coefficient
Gamma0 = (1-beta*H0)./(1+beta*H0);
Gamma = resonatorprofile(nu,nu0,Qu,beta,'voltagereflection');

ok(2) = areequal(Gamma0,Gamma,1e-12,'abs');

% Voltage reflection coefficient
Prefl0 = abs(Gamma).^2;
Prefl = resonatorprofile(nu,nu0,Qu,beta,'powerreflection');

ok(3) = areequal(Prefl0,Prefl,1e-12,'abs');



