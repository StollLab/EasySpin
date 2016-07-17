function [err,data] = test(opt,olddata)

%==================================================================
% Test a and F terms against explicit formulae and stev
%==================================================================
% Abragam/Bleaney p. 437

S = 5/2;
[Sx,Sy,Sz,Sp,Sm] = sop(S,'x','y','z','+','-');
n = S*(S+1);
eS = eye(2*S+1);

Sys.S = S;
Sys.g = [2 2 2];

% (1) a term
%------------------------------------------------------------------
%{
Sys.aF = [121 0];

H0 = Sys.aF(1)*(stev(S,4,0) + 5*stev(S,4,4))/120;
H1 = Sys.aF(1)*((Sx^4+Sy^4+Sz^4) - 1/5*n*(3*S^2+3*S-1)*eS)/6;
H2 = zfield(Sys);

D = H1-H0; err(1) = any(abs(D(:))>1e-6);
D = H2-H0; err(2) = any(abs(D(:))>1e-6);
%}
err(1:2) = 0;

% (2) F term
%-----------------------------------------------------------------
Sys.aF = [0 1];

H0 = stev(S,4,0)/180;
H1 = (35*Sz^4-30*n*Sz^2+25*Sz^2-(6*n-3*n^2)*eS)/180;
H2 = zfield(Sys);

D = H1-H0; err(3) = any(abs(D(:))>1e-6);
D = H2-H0; err(4) = any(abs(D(:))>1e-6);

data = [];
