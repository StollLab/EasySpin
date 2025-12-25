function ok = test()

% Compare resonatorprofile() output with mathematical expression

nu = linspace(9.2,9.8,101); % GHz
nu0 = 9.5; % GHz
Qu = 2000;
beta = 1;

Gamma2 = resonatorprofile(nu,nu0,Qu,beta,'powerreflection');

xi = (nu/nu0-nu0./nu)/2;
Gamma2_ref = 1 - 4*beta./((1+beta)^2+(2*Qu*xi).^2);

ok = areequal(Gamma2,Gamma2_ref,1e-10,'abs');
