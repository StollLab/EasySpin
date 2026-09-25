function ok = test()

% Exp.CenterSweep has precedence over Exp.Range

Exp.CenterSweep = [340 80];
Exp.Range = [100 200];
R1 = runprivate('p_sweeprange',Exp,false,true);

Exp = struct('Range',[300 380]);
R2 = runprivate('p_sweeprange',Exp,false,true);

Exp = struct('mwCenterSweep',[9.5 1],'mwRange',[1 2]);
R3 = runprivate('p_sweeprange',Exp,true,false);

Exp = struct('CenterSweep',NaN,'Range',NaN);
R4 = runprivate('p_sweeprange',Exp,false,true);

ok = isequal(R1,[300 380]) && isequal(R2,[300 380]) && ...
  isequal(R3,[9 10]) && isempty(R4);
