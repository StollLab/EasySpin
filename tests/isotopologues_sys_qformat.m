function ok = test()

% Isotopologues with isotope-specific rescaling of quadrupole couplings

% axial Q
Qlist{1} = [1];  % axial
Qlist{2} = [1 0.3]; % rhombic [eeqQ/h eta]
Qlist{3} = [1 2 4]; % principal values
Qlist{4} = [3 1 2; 4 5 6; 7 4 3]; % full tensor
Qlist{5} = [1 2 -3 0.4 0.5 0.6]; % symmetric matrix [xx yy zz xy xz yz]

qmratio = nucqmom('65Cu')/nucqmom('63Cu');
for k = 1:numel(Qlist)
  Q_63Cu = Qlist{k};
  if numel(Q_63Cu)==2
    Q_65Cu = [Q_63Cu(1)*qmratio Q_63Cu(2)];  % eta is not scaled
  else
    Q_65Cu = Q_63Cu*qmratio;
  end
  Sys.Nucs = 'Cu';
  Sys.Q = Q_63Cu;
  Iso = isotopologues(Sys);
  ok(k) = ...
    areequal(Iso(1).Q,Q_63Cu,1e-10,'rel') && ...
    areequal(Iso(2).Q,Q_65Cu,1e-10,'rel');
end
