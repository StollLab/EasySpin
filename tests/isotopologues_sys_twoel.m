function ok = test()

% Test whether isotopologue() handles multi-electron system

Sys.S = [1/2 1];
Sys.g = [2 2.2];
Sys.Nucs = 'H,H';
Sys.A = [10 10 20, 0 0 0; 0 0 0, 8 8 3];

Iso = isotopologues(Sys);

ok = true;
for k = 1:numel(Iso)
  ok = ok && all(size(Iso(k).A)==[2 6]);
end
