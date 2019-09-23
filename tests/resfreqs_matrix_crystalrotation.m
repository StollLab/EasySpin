function [err,data] = test(opt,olddata)

% Check whether transition intensity follows crystal rotation

clear Sys Exp

Sys.S = 1;
Sys.D = 30e3*2*[1 0.1];
Opt.Transitions = [1 2; 1 3; 2 3];

% Rotation of crystal around z, compensated by B1 rotation around z
for k=1:10
  Exp.CrystalOrientation = [0 0 0];
  Exp.Mode = [0 0];
  [dum,i0] = resfreqs_matrix(Sys,Exp,Opt);
  
  a = rand*pi;
  Exp.CrystalOrientation = [a 0 0];
  Exp.Mode = [0 -a];
  [dum,i1] = resfreqs_matrix(Sys,Exp,Opt);
  
  err1(k) = ~areequal(i0,i1,1e-8,'abs');
end

err = any(err1);
data = [];
