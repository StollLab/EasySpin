function ok = test(opt)

%-------------------------------------------------------------
% ENDOR spectrum of a S=3/2, I=1/2 system
%-------------------------------------------------------------

Sys.S = 3/2;
Sys.D = 200;
Sys.Nucs='1H';
Sys.A = [3 5];
Sys.lwEndor = 0.2;

Exp.Field = 350;
Exp.mwFreq = 9.5;
Exp.Range = [0 30];

Opt.Method = 'matrix';
[nu,spc_matrix] = salt(Sys,Exp,Opt);
Opt.Method = 'perturb2';
[nu,spc_pt2] = salt(Sys,Exp,Opt);

if opt.Display
  plot(nu,spc_matrix,nu,spc_pt2);
  legend('matrix','perturb2');
end

ok = areequal(spc_pt2,spc_matrix,0.01,'rel');
