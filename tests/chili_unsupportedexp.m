function ok = test()

% Unsupported Exp fields are rejected

Sys.g = [2 2.1 2.2];
Sys.lw = 1;
Sys.tcorr = 1e-9;
Exp0.mwFreq = 9.5;
Exp0.Range = [300 350];

fields = {'lightBeam','MolFrame','CrystalSymmetry','SampleRotation','SampleFrame'};
values = {'perpendicular',[0 pi/4 0],'C2h',{'x',pi/4},[0 0 0; 0 pi/4 0]};

for k = numel(fields):-1:1
  Exp = Exp0;
  Exp.(fields{k}) = values{k};
  try
    chili(Sys,Exp);
    ok(k) = false;
  catch err
    ok(k) = contains(err.message,['Exp.' fields{k}]);
  end
end

end
