function [err,data] = test(opt,olddata)

%======================================================
% Syntax test, 2 nuclei
%======================================================
Sys.Nucs = '1H,13C';
Sys.A = [10 3];

% Isotropic coupling
Sys.nn = 0.1;

nnint(Sys);
nnint(Sys,[1 2]);
H = nnint(Sys);
H = nnint(Sys,[1 2]);

% Principal values
Sys.nn = 0.1*[1 1 1];

nnint(Sys);
nnint(Sys,[1 2]);
H = nnint(Sys);
H = nnint(Sys,[1 2]);

% Full tensor
Sys.nn = 0.1*eye(3);

nnint(Sys);
nnint(Sys,[1 2]);
H = nnint(Sys);
H = nnint(Sys,[1 2]);

err = false;
data = [];
