function ok = test()

% Crystal with several sites and a nucleus: output arrays have one row
% per transition and site

Sys.g = [2 2.1 2.2];
Sys.Nucs = '1H';
Sys.A = [10 20 30];  % MHz
Sys.HStrain = [5 6 7];  % MHz

Exp.Field = 350;
Exp.SampleFrame = [0 0 0; 0.1 0.2 0.3];
Exp.MolFrame = [0.4 0.5 0.6];
Exp.CrystalSymmetry = 'P222';  % 4 sites

[Pos,Int,Wid] = resfreqs_perturb(Sys,Exp);
nRows = 2*4;  % 2 nuclear sublevels, 4 sites

ok(1) = isequal(size(Pos),[nRows 2]);
ok(2) = isequal(size(Int),size(Pos));
ok(3) = isequal(size(Wid),size(Pos));

% compare with single-site calculations
for iSite = 4:-1:1
  Opt.Sites = iSite;
  Pos1(:,iSite,:) = resfreqs_perturb(Sys,Exp,Opt);
end
ok(4) = areequal(Pos,reshape(Pos1,nRows,2),1e-10,'rel');
