function ok = test()

% Compare spin systems read from main output files of identical
% calculations run with different ORCA versions. Results differ slightly
% between versions (numerics, spin-orbit treatment), so tolerances are
% loose, but tight enough to catch errors in units, signs, and frames.

ver1 = 'v4.2.1';
ver2 = 'v5.0.4';
mols = {'aminoxyl','dioxygen','nitroxide','tripletformaldehyde'};

% Tolerances
tol.g = 1e-3;
tol.D = 0.1;  % relative to largest principal value
tol.A = 2;  % MHz
tol.Q = 0.01;  % MHz

T = @(v,ang) erot(ang).'*diag(v)*erot(ang);
maxdiff = @(X,Y) max(abs(X(:)-Y(:)));

for m = numel(mols):-1:1
  Sys1 = orca2easyspin(['orca/' ver1 '/' mols{m} '.out']);
  Sys2 = orca2easyspin(['orca/' ver2 '/' mols{m} '.out']);

  okk = Sys1.S==Sys2.S && isequal(Sys1.xyz,Sys2.xyz);
  okk = okk && maxdiff(T(Sys1.g,Sys1.gFrame),T(Sys2.g,Sys2.gFrame))<tol.g;
  okk = okk && isfield(Sys1,'D')==isfield(Sys2,'D');
  if isfield(Sys1,'D')
    Dmax = max(abs(Sys1.D));
    okk = okk && maxdiff(T(Sys1.D,Sys1.DFrame),T(Sys2.D,Sys2.DFrame))<tol.D*Dmax;
  end
  okk = okk && strcmp(Sys1.Nucs,Sys2.Nucs) && isequal(Sys1.NucsIdx,Sys2.NucsIdx);
  for n = 1:numel(Sys1.NucsIdx)
    okk = okk && maxdiff(T(Sys1.A(n,:),Sys1.AFrame(n,:)),T(Sys2.A(n,:),Sys2.AFrame(n,:)))<tol.A;
    okk = okk && maxdiff(T(Sys1.Q(n,:),Sys1.QFrame(n,:)),T(Sys2.Q(n,:),Sys2.QFrame(n,:)))<tol.Q;
  end
  ok(m) = okk;
end
