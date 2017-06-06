function [err,data] = test(opt,olddata)
N =4;
Sys.S = 1/2*ones(1,N);
pairs = nchoosek(N,2);
Sys.ee = rand(pairs,1);
Sys.ee2 = rand(pairs,1);

Sys2.S = [Sys.S,2*Sys.S];
comb1 = nchoosek(1:N,2);
comb2 = nchoosek(1:(2*N),2);
pairs2 = nchoosek(2*N,2);
Sys2.ee = zeros(pairs2,1);
Sys2.ee2 = zeros(pairs2,1);

for n =1:pairs
  ind = logical((comb2(:,1) == comb1(n,1)) .* (comb2(:,2) == comb1(n,2)));
  Sys2.ee(ind) = Sys.ee(n);
  Sys2.ee2(ind) = Sys.ee2(n);
end
% all aditional couplings in the bigger Sys2 are zero, unique eigenvalues 
% have to be identical 
E1 = eig(eeint(Sys));
E2 = eig(eeint(Sys2,1:N));
E1 = unique(round(E1*1e9)/1e9);
E2 = unique(round(E2*1e9)/1e9);

ind2 = (Sys2.ee == 0);
li = sum(ind2);
Sys2.ee(ind2) = rand(li,1);
Sys2.ee2(ind2) = rand(li,1);

% additional couplings are non-zero, but we consider only the smaller
% subsystem, therefore unique eigenvalues have to be identical
E3 = eig(eeint(Sys));
E4 = eig(eeint(Sys2,1:N));
E3 = unique(round(E3*1e9)/1e9);
E4 = unique(round(E4*1e9)/1e9);

err = ~all([areequal(E1,E2),areequal(E3,E4)]);
data = [];