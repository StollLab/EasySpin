% Quick way to compute spin matrices by hand
%==========================================================================
% Matrices are constructed using diag().

% Spin quantum number, change as needed
J = 3/2;

% The simplest one: Jz is diagonal
Jz = diag(J:-1:-J);

% Up- and downshift matrices J+, J-: non-zero
% on upper and lower diagonal, respectively.
values = sqrt((1:2*J).*(2*J:-1:1));
Jp = diag(values,+1)
Jm = diag(values,-1)

% Cartesian spin matrices Jx and Jy are
% linear combinations of J+ and J-
Jx = (Jp + Jm)/2
Jy = (Jp - Jm)/2i
