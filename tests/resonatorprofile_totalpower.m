function ok = test()

% Make sure that reflected plus transmitted power equals 1.

nu = linspace(9.2,9.8,101); % GHz
nu0 = 9.5; % GHz
Qu = 200;
beta = 1;

P_refl = resonatorprofile(nu,nu0,Qu,beta,'powerreflection');
P_trans = resonatorprofile(nu,nu0,Qu,beta,'powertransmission');
P_total = P_refl + P_trans;

ok = areequal(P_total,ones(size(P_total)),1e-15,'abs');
