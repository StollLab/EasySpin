function [err,data] = test(opt,olddata)


% full D matrices with Euler angles
%================================================

Sys.S = 3/2;

D = rand(3);
DFrame = rand(1,3)*2*pi;
Sys.D = D;
Sys.DFrame = DFrame;

H1 = zfield(Sys);

R = erot(DFrame);
Dr = R.'*D*R;
Sys.D = Dr;
Sys.DFrame = [0 0 0];
H2 = zfield(Sys);

deviation = max(abs(H1(:)-H2(:)))/max(abs(H1(:)));


err = deviation>1e-7;
data = [];
