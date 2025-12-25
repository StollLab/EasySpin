function ok = test()

% Compare resonatorprofile() outputs with the mathematical expressions
% that are intended to be implemented in the function.

% Settings
nu = linspace(9,10,1001); % GHz
nu0 = 9.5; % GHz
Qu = 200;
beta = 1.3;

QL = Qu/(1+beta);
xi = nu/nu0-nu0./nu;

% Reference quantities
H = 1./(1+1i*QL*xi);
Gamma = (1-beta+1i*Qu*xi)./(1+beta+1i*Qu*xi);
Gamma_squared = abs(Gamma).^2;
S11 = Gamma;
RL = -20*log10(abs(Gamma));
VSWR = (1+abs(Gamma))./(1-abs(Gamma));

% Transfer function
H_ = resonatorprofile(nu,nu0,Qu,beta,'transferfunction');
ok(1) = areequal(H_,H,1e-12,'abs');

% Voltage reflection coefficient
Gamma_ = resonatorprofile(nu,nu0,Qu,beta,'voltagereflection');
ok(2) = areequal(Gamma_,Gamma,1e-12,'abs');

% Power reflection coefficient
Gamma2_ = resonatorprofile(nu,nu0,Qu,beta,'powerreflection');
ok(3) = areequal(Gamma2_,Gamma_squared,1e-12,'abs');

% S11
S11_ = resonatorprofile(nu,nu0,Qu,beta,'S11');
ok(4) = areequal(S11_,S11,1e-12,'abs');

% Return loss
RL_ = resonatorprofile(nu,nu0,Qu,beta,'returnloss');
ok(5) = areequal(RL_,RL,1e-12,'abs');

% VSWR
VSWR_ = resonatorprofile(nu,nu0,Qu,beta,'VSWR');
ok(6) = areequal(VSWR_,VSWR,1e-12,'abs');
