function ok = test()

DefRelTime = 1e10;

% First Test ---------------------------------------

Sys.Spins = [1/2 1/2];
Sys.T1 = 1.5;
Sys.T2 = 0.5;

RefGamma1 = zeros(16,16);

for i = [2:5, 7:10, 12:15]
  RefGamma1(i,i) = 1/Sys.T2;
end


for i = [1 6 11 16]
  for ii = [1 6 11 16]
    if i ~= ii
      RefGamma1(ii,i) = -1/Sys.T1;
      RefGamma1(i,ii) = -1/Sys.T1;
    end
  end
end


[Gamma1] = runprivate('s_relaxationsuperoperator',Sys);

% Second Test - Default Values for T1 --------------

RefGamma2 = RefGamma1;

for i = [1 6 11 16]
  for ii = [1 6 11 16]
    if i ~= ii
      RefGamma2(ii,i) = -1/DefRelTime;
      RefGamma2(i,ii) = -1/DefRelTime;
    end
  end
end

Sys.T1 = 0; % this value is being set in s_propagationsetup if the field is empty or not found

[Gamma2] = runprivate('s_relaxationsuperoperator',Sys);


% Third Test - Default Values for T2 --------------
Sys.T1 = 1.5;
Sys.T2 = 0;

RefGamma3 = RefGamma1;

for i = [2:5, 7:10, 12:15]
  RefGamma3(i,i) = 1/DefRelTime;
end

[Gamma3] = runprivate('s_relaxationsuperoperator',Sys);

ok = all([isequal(RefGamma1,Gamma1) isequal(RefGamma2,Gamma2) isequal(RefGamma3,Gamma3)]);

