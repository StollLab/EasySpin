function ok = test()

% Compare plegendre results to explicit expressions for low orders L
% Explicit expressions are from Mathematica, LegendreP[L,M,z]

z = 2*rand(1,100)-1;

LM{ 1} = [0, 0]; P{ 1} = ones(size(z));

LM{ 2} = [1,-1]; P{ 2} = sqrt(1-z.^2)/2;
LM{ 3} = [1, 0]; P{ 3} = z;
LM{ 4} = [1, 1]; P{ 4} = -sqrt(1-z.^2);

LM{ 5} = [2,-2]; P{ 5} = -1/8*(z-1).*(z+1);
LM{ 6} = [2,-1]; P{ 6} = 1/2*z.*sqrt(1-z.^2);
LM{ 7} = [2, 0]; P{ 7} = (3*z.^2-1)/2;
LM{ 8} = [2, 1]; P{ 8} = -3*z.*sqrt(1-z.^2);
LM{ 9} = [2, 2]; P{ 9} = -3*(z-1).*(z+1);

LM{10} = [3,-3]; P{10} = (1-z.^2).^(3/2)/48;
LM{11} = [3,-2]; P{11} = -z.*(z.^2-1)/8;
LM{12} = [3,-1]; P{12} = sqrt(1-z.^2).*(5*z.^2-1)/8;
LM{13} = [3, 0]; P{13} = (5*z.^3-3*z)/2;
LM{14} = [3, 1]; P{14} = -3/2*sqrt(1-z.^2).*(5*z.^2-1);
LM{15} = [3, 2]; P{15} = -15*z.*(z.^2-1);
LM{16} = [3, 3]; P{16} = -15*(1-z.^2).^(3/2);

thr = 1e-10;

for k = 1:numel(P)
  L_ = LM{k}(1);
  M_ = LM{k}(2);
  P_ = plegendre(L_,M_,z,true);
  ok(k) = areequal(P_,P{k},thr,'abs');
end
