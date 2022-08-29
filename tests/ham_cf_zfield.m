function ok = test()

% Orbital angular momenta can also be defined as spins, therefore the two
% Hamiltonians should be identical.

SysSL.S = [1/2 1];
SysSL.L = [2 1];
SysSL.soc = zeros(2,1);
SysSL.ee = zeros(nchoosek(length(SysSL.S),2),1);

SysS.S = [SysSL.S,SysSL.L];
SysS.ee = zeros(nchoosek(length(SysS.S),2),1);
for k = 2:2:8
  fieldname = sprintf('CF%d',k);
  spinfieldname = sprintf('B%d',k);
  SysSL.(fieldname) = rand(2,2*k+1).*repmat(((k/2)<=SysSL.L).',1,2*k+1);
  SysS.(spinfieldname) = [zeros(2,2*k+1);SysSL.(fieldname)]; 
end

H_S = ham_zf(SysS);
H_SL = ham_cf(SysSL);

ok = areequal(H_S,H_SL,1e-10,'abs');
