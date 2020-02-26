function ok = test()

% Check number of isotopologues, their abundances, A conversion factors
% and Q conversion factors for a system with 12 equivalent Cl with natural
% abundance

n = 12;
isoList = isotopologues('Cl',n,[],0);

ok = numel(isoList)==13;

% Natural abundances of chlorine isotopes
a1 = nucabund('35Cl');
a2 = nucabund('37Cl');

n2 = 0:n;
n1 = n-n2;

% Calculate abundances
for k = 1:numel(isoList)
  w0(k) = a1^n1(k)*a2^n2(k);
  w0(k) = w0(k)*factorial(n1(k)+n2(k))/factorial(n1(k))/factorial(n2(k));
end

ok = ok && areequal([isoList.weight],w0,1e-10,'rel');

