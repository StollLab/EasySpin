function ok = test()

% D strain (correlated and uncorrelated) with hyperfine coupling, compared
% with resfreqs_matrix

Sys.S = 1;
Sys.g = [2 2.01 2.02];
Sys.D = [300 50];  % MHz
Sys.Nucs = '14N';
Sys.A = [20 25 60];  % MHz
Sys.DStrain = [30 10];  % MHz

Exp.Field = 3350;
Exp.SampleFrame = [0.3 0.7 0.2; 1 0.4 2];

Opt.Threshold = 1e-3;
rDE = [0 0.5];
for k = 1:numel(rDE)
  Sys.DStrainCorr = rDE(k);
  [Pos_m,Int_m,Wid_m] = resfreqs_matrix(Sys,Exp,Opt);
  [Pos_p,~,Wid_p] = resfreqs_perturb(Sys,Exp);
  [dPos,dWid] = compare(Pos_m,Int_m,Wid_m,Pos_p,Wid_p);
  ok(k) = dPos<2 && dWid<0.03;  % MHz, relative
end

end

function [dPos,dWid] = compare(Pos_m,Int_m,Wid_m,Pos_p,Wid_p)
% match each perturbation line to the closest allowed matrix line
dPos = 0;
dWid = 0;
for iOri = 1:size(Pos_p,2)
  allowed = Int_m(:,iOri)>0.3*max(Int_m(:,iOri));
  P = Pos_m(allowed,iOri);
  W = Wid_m(allowed,iOri);
  for r = 1:size(Pos_p,1)
    [d,j] = min(abs(P-Pos_p(r,iOri)));
    dPos = max(dPos,d);
    dWid = max(dWid,abs(W(j)-Wid_p(r,iOri))/W(j));
  end
end
end
