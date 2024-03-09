function ok = test()

% Syntax test
nu = linspace(9,10,1001);
nu0 = 9.6;
Qu = 2000;
beta = 2;

H = resonatorprofile(nu,nu0,Qu,beta,'transferfunction');
H = resonatorprofile(nu,nu0,Qu,'transferfunction');
G = resonatorprofile(nu,nu0,Qu,beta,'voltagereflection');
G = resonatorprofile([],nu0,Qu,beta,'voltagereflection');
P = resonatorprofile(nu,nu0,Qu,beta,'powerreflection');

ok = true;

