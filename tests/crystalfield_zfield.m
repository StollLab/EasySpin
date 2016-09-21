function [err,data] = test(opt,olddata)
%orbital angular momenta can also be defined as spins, therefore the two
%Hamiltonians should be identical

n = randi(3);
Sys.S = randi(3,1,n)/2;
Sys.L = randi(3,1,n);
Sys.soc = zeros(n,1);
if n>1, Sys.ee = zeros(nchoosek(length(Sys.S),2),1);end
PureSpin.S = [Sys.S,Sys.L];
PureSpin.ee = zeros(nchoosek(length([Sys.S,Sys.L]),2),1);
for k=2:2:8
  fieldname = sprintf('CF%d',k);
  spinfieldname = sprintf('B%d',k);
  Sys.(fieldname) = rand(n,2*k+1).*repmat(((k/2)<=Sys.L).',1,2*k+1);
  PureSpin.(spinfieldname) = [zeros(n,2*k+1);Sys.(fieldname)]; 
end

H1 = zfield(PureSpin);
H2 = crystalfield(Sys);

err = ~areequal(H1,H2,1e-10);
data = [];

