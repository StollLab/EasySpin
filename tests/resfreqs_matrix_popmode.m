function ok = test()

% Test PopMode = 'highfield'

clear
Sys.S = 1;
Sys.D = -3*30e3*[1 0.000001];

Exp.Temperature = 0;

[nu1, int1] = resfreqs_matrix(Sys,Exp);

Sys.PopMode = 'highfield';
Sys.Pop = [1 0 0];

[nu2, int2] = resfreqs_matrix(Sys,Exp);

ok = areequal(int1,int2,1e-9,'abs');
