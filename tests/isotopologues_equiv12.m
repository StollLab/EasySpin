function [err,data] = test(opt,olddata)

% Check number of isotopologues, their abundances, A conversion factors
% and Q conversion factors for a system with 12 equivalent Cl with natural
% abundance
%--------------------------------------------------------------------------

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

ok = ok && areequal([isoList.weight],w0);

%{
% Calculate A and Q scaling factor
gn1 = nucgval('35Cl');
gn2 = nucgval('37Cl');
Q1 = nucqmom('35Cl');
Q2 = nucqmom('37Cl');
rA = gn2/gn1;
rQ = Q2/Q1;
Ascale{1} = 1;
Ascale{13} = rA;
Qscale{1} = 1;
Qscale{13} = rQ;
for k = 2:12
  Ascale{k} = [1 rA];
  Qscale{k} = [1 rQ];
end

for k = 1:numel(isoList)
  ok = ok && all(isoList(k).Ascale==Ascale{k});
  ok = ok && all(isoList(k).Qscale==Qscale{k});
end
%}

err = ~ok;
data = [];
