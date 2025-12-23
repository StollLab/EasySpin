function ok = test()

% Syntax test

% Settings
nu = linspace(9,10,1001);  % GHz
nu0 = 9.6;  % GHz
Qu = 2000;
beta = 2;

% Missing frequency axis
H = resonatorprofile([],nu0,Qu,beta,'transferfunction');

% Various output quantities
H = resonatorprofile(nu,nu0,Qu,beta,'transferfunction');
G = resonatorprofile(nu,nu0,Qu,beta,'voltagereflection');
G2 = resonatorprofile(nu,nu0,Qu,beta,'powerreflection');
S11 = resonatorprofile(nu,nu0,Qu,beta,'S11');
RL = resonatorprofile(nu,nu0,Qu,beta,'returnloss');
VSWR = resonatorprofile(nu,nu0,Qu,beta,'VSWR');

ok = true;
