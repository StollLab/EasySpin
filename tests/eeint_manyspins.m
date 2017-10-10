function [err,data] = test(opt,olddata)

%======================================================
% Test for correct ordering of spin pairs
%======================================================

% Generate spin system with 3 coupled spins
N = 3;
Sys1.S = 1/2*ones(1,N);
pairs = nchoosek(N,2);
Sys1.ee = rand(pairs,1);
Sys1.ee2 = rand(pairs,1);

% Add 3 more spins, without coupling to the first 3
Sys2.S = [Sys1.S,Sys1.S];
comb1 = nchoosek(1:N,2);
comb2 = nchoosek(1:(2*N),2);
pairs2 = nchoosek(2*N,2);
Sys2.ee = zeros(pairs2,1);
Sys2.ee2 = zeros(pairs2,1);

% Copy interaction parameters
for n =1:pairs
  ind = logical((comb2(:,1) == comb1(n,1)) .* (comb2(:,2) == comb1(n,2)));
  Sys2.ee(ind) = Sys1.ee(n);
  Sys2.ee2(ind) = Sys1.ee2(n);
end

% Generat e-e interaction Hamiltonians
H1 = eeint(Sys1);
H2 = eeint(Sys2);

% Determine whether there are elements that are different
valdiff = setdiff(H1(:),H2(:));

err = ~isempty(valdiff);
data = [];
