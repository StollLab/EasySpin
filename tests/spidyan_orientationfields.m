function ok = test()

% spidyan simulates a single orientation and rejects sample orientation fields

Sys.S = 1/2;
Sys.ZeemanFreq = 9.5;  % GHz

Pulse.Type = 'rectangular';
Pulse.tp = 0.016;  % µs
Pulse.Flip = pi/2;  % rad

Exp.Sequence = {Pulse 0.1};
Exp.mwFreq = 9.5;  % GHz

fields = {'SampleFrame','MolFrame','CrystalSymmetry','SampleRotation'};
values = {[0 1 0],[0 1 0],'C2h',{1,pi/4}};
for k = 1:numel(fields)
  Exp_ = Exp;
  Exp_.(fields{k}) = values{k};
  try
    spidyan(Sys,Exp_);
    ok(k) = false;
  catch err
    ok(k) = contains(err.message,fields{k});
  end
end
