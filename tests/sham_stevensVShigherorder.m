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

Sys1.S = randi(10)/2+4;
% l = lB+lS = lS-lB, m =l, ..., -l
Sys1.ZB02.l = 2;
Sys1.ZB02.vals = rand(1,5);

Sys1.ZB04.l=4;
Sys1.ZB04.vals = rand(1,9);

Sys1.ZB06.l=6;
Sys1.ZB06.vals = rand(1,13);

Sys1.ZB08.l=8;
Sys1.ZB08.vals = rand(1,17);

Sys2.S = Sys1.S;
%conversion to Stevens opertor parameters
Sys2.B2 = Sys1.ZB02.vals./Alms(2,1:5);
Sys2.B4 = Sys1.ZB04.vals./Alms(4,1:9);
Sys2.B6 = Sys1.ZB06.vals./Alms(6,1:13);
Sys2.B8 = Sys1.ZB08.vals./Alms(8,1:17);

%obtain Hamiltonian matrices
Hz = sham(Sys1,[0,0,0]);
H = sham(Sys2,[0,0,0]);

% test
threshold = 1e-6;
err = ~areequal(H,Hz,threshold);
data =[];



