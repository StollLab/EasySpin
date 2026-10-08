function ok = test()

% Test strains in resfreqs_matrix 

% HStrain
%-------------------------------------------------------------------------------
clear Sys Exp
Exp.SampleFrame = rand(1,3)*2*pi;
R = erot(Exp.SampleFrame);
Sys.HStrain = [1 0 0];
[p,i,w] = resfreqs_matrix(Sys,Exp);
ok(1) = areequal(w,abs(R(1,3)),1e-5,'abs');

clear Sys Exp
Exp.SampleFrame = rand(1,3)*2*pi;
R = erot(Exp.SampleFrame);
Sys.HStrain = [0 1 0];
[p,i,w] = resfreqs_matrix(Sys,Exp);
ok(2) = areequal(w,abs(R(2,3)),1e-5,'abs');

clear Sys Exp
Exp.SampleFrame = rand(1,3)*2*pi;
R = erot(Exp.SampleFrame);
Sys.HStrain = [0 0 1];
[p,i,w] = resfreqs_matrix(Sys,Exp);
ok(3) = areequal(w,abs(R(3,3)),1e-5,'abs');

% D strain
%-------------------------------------------------------------------------------
clear Sys Exp
Sys.S = 3/2;
Sys.D = rand*1000;
DS = rand*Sys.D;
Sys.StrainPars = {'D'};
Sys.StrainFWHM = DS;
Exp.Field = rand*1000;
[p,i,w] = resfreqs_matrix(Sys,Exp);
ok(4) = any(abs(w)<1e-8*DS);
ok(5) = any(abs(w-2*DS)<1e-8*DS);




