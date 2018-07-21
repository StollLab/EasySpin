function [err,data] = test(opt,olddata)

DefRelTime = 1e10;

% First Test ---------------------------------------

Sys.Spins = 1/2;
Sys.T1 = [0 1.5; 1 0];
Sys.T2 = 0.5;

RefGamma1(2,2) = 1/Sys.T2;
RefGamma1(3,3) = 1/Sys.T2;
RefGamma1(1,4) = -1/Sys.T1(1,2);
RefGamma1(4,1) = -1/Sys.T1(2,1);

[Gamma1] = runprivate('s_relaxationsuperoperator',Sys);

% Second Test - Have Default for T2 --------------

Sys.T2 = 0; % this value is being set in s_propagationsetup if the field is empty or not found

RefGamma2 = RefGamma1;
RefGamma2(2,2) = 1/DefRelTime;
RefGamma2(3,3) = 1/DefRelTime;

[Gamma2] = runprivate('s_relaxationsuperoperator',Sys);

% Third Test - Matrix for T2 --------------
Sys.T1 = [0 1.5; 0 0];
Sys.T2 = [0 0.5; 0 0];

RefGamma3(2,2) = 1/Sys.T2(1,2);
RefGamma3(3,3) = 1/Sys.T2(1,2);
RefGamma3(1,4) = -1/Sys.T1(1,2);
RefGamma3(4,1) = -1/Sys.T1(1,2);

[Gamma3] = runprivate('s_relaxationsuperoperator',Sys);

% Third Test - Empty matrix for T1 --------------
Sys.T1 = [0 0; 0 0];
Sys.T2 = [0 0.5; 0.3 0]; % checks if the value 0.3 is indeed ignored

RefGamma4(2,2) = 1/Sys.T2(1,2);
RefGamma4(3,3) = 1/Sys.T2(1,2);
RefGamma4(1,4) = -1/DefRelTime;
RefGamma4(4,1) = -1/DefRelTime;

[Gamma4] = runprivate('s_relaxationsuperoperator',Sys);

if any([~isequal(RefGamma1,Gamma1) ~isequal(RefGamma2,Gamma2) ~isequal(RefGamma3,Gamma3)  ~isequal(RefGamma4,Gamma4)])
  err = 1;
else
  err = 0;
end

data = [];

