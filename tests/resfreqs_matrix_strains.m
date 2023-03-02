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

% DStrain
%-------------------------------------------------------------------------------
clear Sys Exp
Sys.S = 3/2;
Sys.D = rand*1000;
Sys.DStrain = rand*Sys.D;
Exp.Field = rand*1000;
[p,i,w] = resfreqs_matrix(Sys,Exp);
ok(4) = ~isempty(find(w==0));
ok(5) = ~isempty(find(w==(2*Sys.DStrain)));




