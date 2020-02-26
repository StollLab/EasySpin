function ok = test()

% vector of angular momentum quantum numbers
J = [1/2 3/2 1 1 1 1];

N = hsdim(J);

Nref = prod(2*J+1);

ok = (N==Nref);
