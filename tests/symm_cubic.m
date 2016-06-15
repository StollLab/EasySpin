function [err,data] = test(opt,olddata)

%======================================================
% Cubic high-spin systems
%======================================================
for n=0:1
  
sys.S = 7/2;
sys.g = 2*[1 1 1];
%------------------------------------------------------
B2 = 0;
B4 = 0;
B6 = 10;

sys.B2 = B2;
sys.B4 = [5*B4 0 0 0 B4 0 0 0 0];
sys.B6 = [0 0 -21*B6 0 0 0 B6 0 0 0 0 0 0];

G = symm(sys);
err(1+n) = ~strcmp(G,'Oh');

%------------------------------------------------------
B2 = 0;
B4 = 10;
B6 = 0;
sys.B2 = B2;
sys.B4 = [5*B4 0 0 0 B4 0 0 0 0];
sys.B6 = [0 0 -21*B6 0 0 0 B6 0 0 0 0 0 0];

G = symm(sys);
err(2+n) = ~strcmp(G,'Oh');

%------------------------------------------------------
B2 = 0;
B4 = 10;
B6 = 12.87;
sys.B2 = B2;
sys.B4 = [5*B4 0 0 0 B4 0 0 0 0];
sys.B6 = [0 0 -21*B6 0 0 0 B6 0 0 0 0 0 0];

G = symm(sys);
err(3+n) = ~strcmp(G,'Oh');
sys.ZB11.l= 0;sys.ZB11.vals= 1; %isotropic term, no change in symmetry
end

data = [];
