function [err,data] = test(opt,olddata)
%orbital angular momenta can also be defined as spins, therefore the two
%Hamiltonians should be identical
n = randi_(3);
Sys.S = randi_(3,1,n)/2;
Sys.L = randi_(2,1,n);
Sys.soc = rand(n,2);
Sys.orf = rand(n,1);
lenS = length(Sys.S);


if n>1, Sys.ee = zeros(nchoosek(lenS,2),1);end
PureSpin.S = [Sys.S,Sys.L];

%distribute soc over ee and ee2
len = 2*lenS;
k = nchoosek(1:len,2);
eelen = nchoosek(len,2);
PureSpin.ee = zeros(eelen,1);
PureSpin.ee2 = zeros(eelen,1);
for m=1:lenS
  x = logical((k(:,1)==m).*(k(:,2)==m+lenS));
  PureSpin.ee(x) = Sys.orf(m)*Sys.soc(m,1);
  PureSpin.ee2(x) = Sys.orf(m)*Sys.orf(m)*Sys.soc(m,2);
end

H1 = eeint(PureSpin);
H2 = soint(Sys);

err = ~areequal(H1,H2,1e-10,'abs');
data = [];

