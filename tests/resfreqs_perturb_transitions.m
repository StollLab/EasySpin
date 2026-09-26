function ok = test()

% Level indices in transition list agree with resfreqs_matrix for S>1/2

Sys.S = 3/2;
Sys.g = 2;
Sys.D = 500;  % MHz

Exp.Field = 350;
Exp.SampleFrame = [0 0 0];

[Pos_p,~,~,Tr_p] = resfreqs_perturb(Sys,Exp);
[Pos_m,~,~,Tr_m] = resfreqs_matrix(Sys,Exp);

[~,ia,ib] = intersect(Tr_m,Tr_p,'rows');
ok(1) = numel(ia)==size(Tr_p,1);
ok(2) = areequal(Pos_m(ia),Pos_p(ib),1e-4,'rel');
