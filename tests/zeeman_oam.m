function [err,data] = test(opt,olddata)
%orbital angular momenta can also be defined as spins, therefore the two
%Hamiltonians should be identical

n = randi(3);
Sys.S = randi(3,1,n)/2;
Sys.L = randi(3,1,n);
Sys.soc = zeros(n,1);
Sys.g = rand(3*n,3);
Sys.orf = rand(n, 1);
if n>1, Sys.ee = zeros(nchoosek(length(Sys.S),2),1);end
PureSpin.S = [Sys.S,Sys.L];
PureSpin.ee = zeros(nchoosek(length([Sys.S,Sys.L]),2),1);
PureSpin.g = [Sys.g;zeros(3*n,3)];
for k=1:n
  PureSpin.g(3*(n+k-1)+1:3*(n+k),:) = -diag(Sys.orf(k)*ones(1,3));
end


H1 = zfield(PureSpin);
H2 = crystalfield(Sys);

err = ~areequal(H1,H2,1e-10);
data = [];