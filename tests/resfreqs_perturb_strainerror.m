function ok = test()

% Perturbation theory: strains of g, A and D are supported by resfields_perturb
% and resfreqs_perturb, but not by endorfrq_perturb. Strains of Q and sigma,
% and strains with Opt.ImmediateBinning, give an error.

Exp.mwFreq = 9.5;
Exp.Field = 330;
Exp.Range = [280 360];

% g strain
Sys.g = [2 2.1 2.2];
Sys.Nucs = '1H';
Sys.A = [10 20 30];
Sys.StrainPars = {'g(1)'};
Sys.StrainFWHM = 0.01;

ok(1) = ~raiseserror(@()resfields_perturb(Sys,Exp));
ok(2) = ~raiseserror(@()resfreqs_perturb(Sys,Exp));
ok(3) = raiseserror(@()endorfrq_perturb(Sys,Exp),'StrainPars');

% Q and sigma strains
SysQ = struct('Nucs','14N','A',[10 20 30],'Q',1,'StrainPars',{{'Q'}},'StrainFWHM',0.1);
SysS = struct('Nucs','1H','A',[10 20 30],'sigma',[1 1 1.001],'StrainPars',{{'sigma(3)'}},'StrainFWHM',1e-4);
ok(4) = raiseserror(@()resfields_perturb(SysQ,Exp),'Sys.Q/Sys.sigma');
ok(5) = raiseserror(@()resfreqs_perturb(SysQ,Exp),'Sys.Q/Sys.sigma');
ok(6) = raiseserror(@()resfields_perturb(SysS,Exp),'Sys.Q/Sys.sigma');
ok(7) = raiseserror(@()resfreqs_perturb(SysS,Exp),'Sys.Q/Sys.sigma');

% Opt.ImmediateBinning
Exp.nPoints = 1024;
Exp.AccumWeights = 1;
Opt.ImmediateBinning = true;
ok(8) = raiseserror(@()callwidths(Sys,Exp,Opt),'ImmediateBinning');

end

%-------------------------------------------------------------------------------
function err = raiseserror(fcn,msg)
try
  fcn();
  err = false;
catch ME
  err = nargin<2 || contains(ME.message,msg);
end
end

function callwidths(Sys,Exp,Opt)
[~,~,~] = resfields_perturb(Sys,Exp,Opt);
end
