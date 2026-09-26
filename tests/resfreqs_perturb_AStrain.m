function ok = test()

% Correlated g/A strain for S=1 (checks ordering of widths across
% electron and nuclear transitions), compared with resfreqs_matrix

Sys.S = 1;
Sys.g = 2;
Sys.D = [200 0];  % MHz
Sys.Nucs = '1H';
Sys.A = [30 30 60];  % MHz
Sys.gStrain = [0.001 0.001 0.002];
Sys.AStrain = [5 5 10];  % MHz

Exp.Field = 3350;
Exp.SampleFrame = [0 0 0];  % B0 along molecular z

Opt.Threshold = 1e-3;
[Pos_m,Int_m,Wid_m] = resfreqs_matrix(Sys,Exp,Opt);
[Pos_p,~,Wid_p] = resfreqs_perturb(Sys,Exp);

allowed = Int_m>0.3*max(Int_m);
Pos_m = Pos_m(allowed);
Wid_m = Wid_m(allowed);
for r = numel(Pos_p):-1:1
  [dPos(r),j] = min(abs(Pos_m-Pos_p(r)));
  dWid(r) = abs(Wid_m(j)-Wid_p(r))/Wid_m(j);
end

ok(1) = max(dPos)<0.1;  % MHz
ok(2) = max(dWid)<0.1;
ok(3) = numel(unique(round(Wid_p)))>1;  % widths differ between lines
