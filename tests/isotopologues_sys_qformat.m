function [err,data] = test(opt,olddata)

% Isotopologues with isotope-specific rescaling of quadrupole couplings
%-------------------------------------------------------------------------------

% axial Q
Qlist{1} = [1];  % axial
Qlist{2} = [1 2]; % rhombic
Qlist{3} = [1 2 4]; % principal values
Qlist{4} = [3 1 2; 4 5 6; 7 4 3]; % full tensor

err = false;
for k = 1:numel(Qlist)
  Q_63Cu = Qlist{k};
  Q_65Cu = Q_63Cu*nucqmom('65Cu')/nucqmom('63Cu');
  Sys.Nucs = 'Cu';
  Sys.Q = Q_63Cu;
  Iso = isotopologues(Sys);
  err = err || ~areequal(Iso(1).Q,Q_63Cu) || ~areequal(Iso(2).Q,Q_65Cu);
end

data = [];
