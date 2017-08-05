function [err,data] = test(opt,olddata)

%======================================================
% Syntax check, three nuclear spins
%======================================================
Sys.Nucs = '1H,13C,15N';
Sys.A = [10 20 30];

% Single value
Sys.nn = [1 2 3];
H = nnint(Sys);
H = nnint(Sys,[1 2]);
H = nnint(Sys,[1 3]);
H = nnint(Sys,[2 3]);

% Principal values
Sys.nn = [1 1 1; 2 2 2; 3 3 3];
H = nnint(Sys);
H = nnint(Sys,[1 2]);
H = nnint(Sys,[1 3]);
H = nnint(Sys,[2 3]);

% Full matrices
Sys.nn = rand(9,3);
H = nnint(Sys);
H = nnint(Sys,[1 2]);
H = nnint(Sys,[1 3]);
H = nnint(Sys,[2 3]);

err = false;
data = [];
