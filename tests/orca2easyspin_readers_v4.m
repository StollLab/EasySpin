function ok = test()

% Compare spin systems read from the main output file and from the
% binary property file (.prop) from ORCA 4,
% for the same ORCA calculations

file1 = 'v4.2.1/%s.out';
file2 = 'v4.2.1/%s.prop';
mols = {'aminoxyl','dioxygen','nitroxide','tripletformaldehyde'};

% Tolerances
tol.xyz = 1e-5;  % Angstrom
tol.g = 1e-6;
tol.D = 1;  % MHz
tol.A = 0.05;  % MHz
tol.Q = 1e-3;  % MHz

for m = numel(mols):-1:1
  Sys1 = orca2easyspin(['orca/' sprintf(file1,mols{m})]);
  Sys2 = orca2easyspin(['orca/' sprintf(file2,mols{m})]);
  ok(m) = compare(Sys1,Sys2,tol);
end

end

%-------------------------------------------------------------------------------
function ok = compare(Sys1,Sys2,tol)

T = @(v,ang) erot(ang).'*diag(v)*erot(ang);
maxdiff = @(X,Y) max(abs(X(:)-Y(:)));

ok = Sys1.S==Sys2.S && Sys1.charge==Sys2.charge;
ok = ok && isequal(Sys1.Elements,Sys2.Elements);
ok = ok && maxdiff(Sys1.xyz,Sys2.xyz)<tol.xyz;
ok = ok && maxdiff(T(Sys1.g,Sys1.gFrame),T(Sys2.g,Sys2.gFrame))<tol.g;
ok = ok && isfield(Sys1,'D')==isfield(Sys2,'D');
if isfield(Sys1,'D')
  ok = ok && maxdiff(T(Sys1.D,Sys1.DFrame),T(Sys2.D,Sys2.DFrame))<tol.D;
end
ok = ok && strcmp(Sys1.Nucs,Sys2.Nucs) && isequal(Sys1.NucsIdx,Sys2.NucsIdx);
for n = 1:numel(Sys1.NucsIdx)
  ok = ok && maxdiff(T(Sys1.A(n,:),Sys1.AFrame(n,:)),T(Sys2.A(n,:),Sys2.AFrame(n,:)))<tol.A;
  ok = ok && maxdiff(T(Sys1.Q(n,:),Sys1.QFrame(n,:)),T(Sys2.Q(n,:),Sys2.QFrame(n,:)))<tol.Q;
end

end
