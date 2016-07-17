function [err,data] = test(opt,olddata)
%==================================================================
% Test whether the 0th order in B of zeemanho is identical (Except the
% known prefactors to results using the Stevens operators
%==================================================================

%Table of conversion factors for Stevens operator
Alm(8,:) = [24*sqrt(1430),2*sqrt(1430),4*sqrt(143/7),2*sqrt(78/7), ...
  4*sqrt(130/7), 2*sqrt(10/7), 4*sqrt(15), 2*sqrt(2), 8*sqrt(2)];
Alm(7,1:8) = [4*sqrt(429), 8*sqrt(429/7),4*sqrt(286/7),8*sqrt(143/7),...
  4*sqrt(13/7), 8*sqrt(13/7),4*sqrt(2/7),8];
Alm(6,1:7) = [4*sqrt(231),2*sqrt(11),4*sqrt(22/5),2*sqrt(22/5),...
  4*sqrt(11/3), 2*sqrt(2/3), 4*sqrt(2)];
Alm(5,1:6) = [6*sqrt(14), 2*sqrt(42/5),sqrt(6/5),12/sqrt(5), 2*sqrt(2/5), 4];
Alm(4,1:5) = [2*sqrt(70), sqrt(7), sqrt(14), 1, 2*sqrt(2)];
Alm(3,1:4) = [sqrt(10), 2*sqrt(5/3), sqrt(2/3),2];
Alm(2,1:3) = [sqrt(6), 1/sqrt(2), sqrt(2)];
Alm(1,1:2) = [1,1];

for n=8:-1:1
    Alms(n,1:2*n+1)=[fliplr(Alm(n,2:n+1)), Alm(n,1:n+1)];
end

Sys.S = ceil(rand*10)/2+4;
Sys2.S = Sys.S ;

lB = 0;
for lS=2:2:8
  Hamstr = ['Ham',num2str([lB, lS,lS],'%i%i%i')];
  Bstr = ['B',num2str(lS,'%i')];
  len = 2*lS+1;
  Sys.(Hamstr) = rand(1,len);
  Sys2.(Bstr) = Sys.(Hamstr)./Alms(lS,1:len);
end


%obtain Hamiltonian matrices
Hz = zeemanho(Sys,[0,0,0]);
H = zfield(Sys2);

% test
threshold = 1e-6;
err = ~areequal(H,Hz,threshold);
data =[];



