function [err,data] = test(opt,olddata)

DefRelTime = 1e10;

% First Test ---------------------------------------

Sys.Spins = [1/2 1/2];
Sys.T1 = [0 1.5 1.5 1.5; 0 0 1.5 1.5; 0 0 0 1.5; 0 0 0 0];
Sys.T2 = 0.5;

RefGamma1 = zeros(16,16);

for i = [2:5, 7:10, 12:15]
  RefGamma1(i,i) = 1/Sys.T2;
end


for i = [1 6 11 16]
  for ii = [1 6 11 16]
    if i ~= ii
      RefGamma1(ii,i) = -1/Sys.T1(1,2);
      RefGamma1(i,ii) = -1/Sys.T1(1,2);
    end
  end
end


[Gamma1] = runprivate('s_relaxationsuperoperator',Sys);

% Second Test - Default Values for T1 --------------
Sys.T1 = [0 1.5 1.5 1.5; 1 0 1.5 1.5; 1 1 0 1.5; 1 1 1 0];
Sys.T2 = 0.5;

RefGamma2 = zeros(16,16);

for i = [2:5, 7:10, 12:15]
  RefGamma2(i,i) = 1/Sys.T2;
end


for i = [1 6 11 16]
  for ii = [1 6 11 16]
    if i ~= ii
      RefGamma2(ii,i) = -1/Sys.T1(1,2);
      RefGamma2(i,ii) = -1/Sys.T1(2,1);
    end
  end
end


[Gamma2] = runprivate('s_relaxationsuperoperator',Sys);


% Third Test - Default Values for T2 --------------
Sys.T2 = [0 .5 .5 .5; 0 0 .5 .5; 0 0 0 .5; 0 0 0 0];

RefGamma3 = RefGamma2;

for i = [2:5, 7:10, 12:15]
  RefGamma3(i,i) = 1/Sys.T2(1,2);
end

[Gamma3] = runprivate('s_relaxationsuperoperator',Sys);

if any([~isequal(RefGamma1,Gamma1) ~isequal(RefGamma2,Gamma2) ~isequal(RefGamma3,Gamma3)])
  err = 1;
else
  err = 0;
end

data = [];

