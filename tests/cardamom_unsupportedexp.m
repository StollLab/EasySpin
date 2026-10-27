function ok = test()

% Unsupported Exp fields are rejected

Sys.g = [2 2.1 2.2];
Sys.lw = 1;
Sys.tcorr = 1e-9;
Exp0.mwFreq = 9.5;
Exp0.Range = [300 350];
Par.Model = 'diffusion';
Par.dtSpin = 1e-9;
Par.nSteps = 10;
Par.nTraj = 1;

fields = {'mwMode','lightBeam','Ordering',...
  'SampleFrame','SampleRotation','MolFrame','CrystalSymmetry'};
values = {'parallel','perpendicular',1,...
  [0 pi/4 0],{'x',pi/4},[0 pi/4 0],'C2h'};

for k = numel(fields):-1:1
  Exp = Exp0;
  Exp.(fields{k}) = values{k};
  try
    cardamom(Sys,Exp,Par);
    ok(k) = false;
  catch err
    ok(k) = contains(err.message,['Exp.' fields{k}]);
  end
end

end
