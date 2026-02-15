function [ok,data] = test(opt,olddata)

% Tests if asymmetric A tensor works

Sys.Nucs = '1H';
Sys.g = [2.2 2.1 2];
Sys.A = 2*[-1.2,-0.8,2];
Sys.lwEndor = [0 0.01];
Exp.Range = [10 18];
Exp.Field = 326.5;
Exp.nPoints = 4000;
Exp.mwFreq = 9.7;

Opt.Threshold = 1e-5;
Opt.GridSize = [20 5];
Opt.Intensity = 'on';
Opt.Enhancement = 'off';
Opt.Verbosity = opt.Verbosity;

Opt.GridSymmetry = 'Ci';
Opt.GridFrame = [0 0 0];

% Diagonal A matrix
Sys.A = diag(Sys.A);
[nu,spc1] = salt(Sys,Exp,Opt);

% A matrix with antisymmetric contribution
Sys.A(1,2) = 5;
Sys.A(2,1) = -5;
[nu,spc2] = salt(Sys,Exp,Opt);

if opt.Display
  subplot(3,1,[1 2]);
  plot(nu,spc1,'b',nu,spc2,'r',nu,olddata.spc1,'b.',nu,olddata.spc2,'r.');
  title('C_i symmetry');
  legend('symmetric A','asymmetric A','symmetric A, ref','asymmetric A, ref');
  xlabel('frequency (MHz)');
  if ~isempty(olddata)
    subplot(3,1,3);
    plot(nu,spc1-olddata.spc1,'b.',nu,spc2-olddata.spc2,'r.');
  end
end

data.spc1 = spc1;
data.spc2 = spc2;

if isempty(olddata)
  ok = [];
else
  ok = areequal(olddata.spc1,spc1,1e-10,'rel') && areequal(olddata.spc2,spc2,1e-10,'rel');
end
