function ok = test()

DefRelTime = 1e10;

% First Test ---------------------------------------

Sys.Spins = 1/2;
Sys.T1 = 1.5;
Sys.T2 = 0.5;

RefGamma1(2,2) = 1/Sys.T2;
RefGamma1(3,3) = 1/Sys.T2;
RefGamma1(1,4) = -1/Sys.T1;
RefGamma1(4,1) = -1/Sys.T1;

[Gamma1] = runprivate('s_relaxationsuperoperator',Sys);

% Second Test - Default Values for T1 --------------

RefGamma2 = RefGamma1;
RefGamma2(1,4) = -1/DefRelTime;
RefGamma2(4,1) = -1/DefRelTime;

Sys.T1 = 0; % this value is being set in s_propagationsetup if the field is empty or not found

[Gamma2] = runprivate('s_relaxationsuperoperator',Sys);


% Third Test - Default Values for T2 --------------
Sys.T1 = 1.5;
Sys.T2 = 0;

RefGamma3 = RefGamma1;
RefGamma3(2,2) = 1/DefRelTime;
RefGamma3(3,3) = 1/DefRelTime;

[Gamma3] = runprivate('s_relaxationsuperoperator',Sys);

ok = all([isequal(RefGamma1,Gamma1) isequal(RefGamma2,Gamma2) isequal(RefGamma3,Gamma3)]);

