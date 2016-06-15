function [err,data] = test(opt,olddata)


% full g matrices with Euler angles
%================================================

B = rand(1,3)*340;

Sys.S = 3/2;
Sys.Nucs = '1H';
Sys.A = [30 40 50];

g = rand(3);
gFrame = rand(1,3)*2*pi;
Sys.g = g;
Sys.gFrame = gFrame;

H1 = zeeman(Sys,B);

R = erot(gFrame);
gr = R.'*g*R;
Sys.g = gr;
Sys.gFrame = [0 0 0];
H2 = zeeman(Sys,B);

deviation = max(abs(H1(:)-H2(:)))/max(abs(H1(:)));


err = deviation>1e-7;
data = [];
