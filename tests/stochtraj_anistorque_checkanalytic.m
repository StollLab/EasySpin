function [err,data] = test(opt,olddata)
% Check that stochtraj_calcanistorque agrees with analytical D-matrix elements from 
% tables for L = 1,2 in Sezer, et al., J.Chem.Phys. 128, 165106 (2008)
% or calculated in Mathematica

rng(1)

tol = 5*eps;
len = 3;

beta  =   pi*(2*rand(1,len)-1);
alpha = 2*pi*(2*rand(1,len)-1);
gamma = 2*pi*(2*rand(1,len)-1);

q = euler2quat(alpha,beta,gamma);

A = q(1,:) - 1i*q(4,:);
B = -q(3,:) - 1i*q(2,:);
Ast = q(1,:) + 1i*q(4,:);
Bst = -q(3,:) + 1i*q(2,:);

ReA = real(A);
ImA = imag(A);
ReB = real(B);
ImB = imag(B);

X = A.*Bst + Ast.*B;
Y = -1i*(A.*Bst - Ast.*B);
Z = A.*Ast - B.*Bst;

% Note that there are typos in Table 3 of Lynden-Bell, Stone, Mol. Sim. 1989,
% here is the correct matrix, which is actually orthogonal:
%
% R = [ 1/2*(A^2-B^2+Ast^2-Bst^2), -1i/2*(A^2+B^2-Ast^2-Bst^2),      -A*B-Ast*Bst;
%      1i/2*(A^2-B^2-Ast^2+Bst^2),   1/2*(A^2+B^2+Ast^2+Bst^2), -1i*(A*B-Ast*Bst);
%                     A*Bst+Ast*B,           -1i*(A*Bst-Ast*B),      A*Ast-B*Bst];

% Test lambda^1_0,0
% -------------------------------------------------------------------------
L = 1;
M = 0;
K = 0;
LMK = [L,M,K];
lambda = 2*rand()-1 + 1i*(2*rand()-1);

numtorque = [ -real(lambda)*Y;
               real(lambda)*X;
               zeros(1,len)];

errl100val = abs(...
                 runprivate('stochtraj_calcanistorque',LMK,[real(lambda),imag(lambda)],q) ...
                 -numtorque...
                );
errl100 = errl100val>tol;

% Test lambda^1_0,1
% -------------------------------------------------------------------------
L = 1;
M = 0;
K = 1;
LMK = [L,M,K];
lambda = 2*rand()-1 + 1i*(2*rand()-1);

numtorque = [                  1/sqrt(2)*imag(lambda)*Z;
                               1/sqrt(2)*real(lambda)*Z;
              -1/sqrt(2)*(real(lambda)*Y+imag(lambda)*X) ];

errl10p1val = abs(...
                  runprivate('stochtraj_calcanistorque',LMK,[real(lambda),imag(lambda)],q)-...
                  numtorque...
                 );
errl10p1 = errl10p1val>tol;

% Test lambda^1_1,0
% -------------------------------------------------------------------------
L = 1;
M = 1;
K = 0;
LMK = [L,M,K];
coef = 1;
lambda = 2*rand()-1 + 1i*(2*rand()-1);

numtorque = [ 1/sqrt(2)*( imag(lambda)*(ReA.^2+ReB.^2-ImA.^2-ImB.^2)+2*real(lambda)*(ReA.*ImA+ReB.*ImB) );
              1/sqrt(2)*( 2*imag(lambda)*(ReA.*ImA-ReB.*ImB) - real(lambda)*(ReA.^2-ReB.^2-ImA.^2+ImB.^2) );
                                                                  zeros(1,len)];

errl1p10val = abs(...
                  runprivate('stochtraj_calcanistorque',LMK,[real(lambda),imag(lambda)],q)-...
                  numtorque...
                  );

errl1p10 = errl1p10val>tol;

% Test lambda^1_1,1
% -------------------------------------------------------------------------
L = 1;
M = 1;
K = 1;
LMK = [L,M,K];
lambda = 2*rand()-1 + 1i*(2*rand()-1);

numtorque = [  real(lambda)*(ReA.*ImB+ImA.*ReB) + imag(lambda)*(ReA.*ReB-ImA.*ImB);
               real(lambda)*(ReA.*ReB-ImA.*ImB) - imag(lambda)*(ReA.*ImB+ImA.*ReB);
               2*real(lambda)*ReA.*ImA  +     imag(lambda)*(ReA.^2-ImA.^2)];

errl1p1p1val = abs(...
                   runprivate('stochtraj_calcanistorque',LMK,[real(lambda),imag(lambda)],q)-...
                   numtorque...
                  );

errl1p1p1 = errl1p1p1val>tol;

% Test lambda^1_-1,1
% -------------------------------------------------------------------------
L = 1;
M = -1;
K = 1;
LMK = [L,M,K];
lambda = 2*rand()-1 + 1i*(2*rand()-1);

numtorque = [  real(lambda)*(ReA.*ImB+ImA.*ReB) - imag(lambda)*(ReA.*ReB-ImA.*ImB);
              -real(lambda)*(ReA.*ReB-ImA.*ImB) - imag(lambda)*(ReA.*ImB+ImA.*ReB);
              -2*real(lambda)*ReB.*ImB +     imag(lambda)*(ReB.^2-ImB.^2) ];

errl1m1p1val = abs(...
                   runprivate('stochtraj_calcanistorque',LMK,[real(lambda),imag(lambda)],q)-...
                   numtorque...
                  );

errl1m1p1 = errl1m1p1val>tol;

% Test lambda^2_0,0
% -------------------------------------------------------------------------
L = 2;
M = 0;
K = 0;
LMK = [L,M,K];
lambda = 2*rand()-1 + 1i*(2*rand()-1);

numtorque = [ 6*real(lambda)*(ReA.*ImB-ImA.*ReB).*Z;
              6*real(lambda)*(ReA.*ReB+ImA.*ImB).*Z;
                                        zeros(1,len) ];

errl200val = abs(...
                 runprivate('stochtraj_calcanistorque',LMK,[real(lambda),imag(lambda)],q)-...
                 numtorque...
                );
errl200 = errl200val>tol;

% Test lambda^2_0,1
% -------------------------------------------------------------------------
L = 2;
M = 0;
K = 1;
LMK = [L,M,K];
lambda = 2*rand()-1 + 1i*(2*rand()-1);

numtorque = [ sqrt(3/2)/2*(2*real(lambda)*X.*Y+imag(lambda)*(3*Z.^2+X.^2-Y.^2-1));
              sqrt(3/2)/2*(2*imag(lambda)*X.*Y+real(lambda)*(3*Z.^2-X.^2+Y.^2-1));
                                    -sqrt(3/2)*(imag(lambda)*X+real(lambda)*Y).*Z ];

errl20p1val = abs(...
                  runprivate('stochtraj_calcanistorque',LMK,[real(lambda),imag(lambda)],q)-...
                  numtorque...
                 );

errl20p1 = errl20p1val>tol;

% Test lambda^2_0,2
% -------------------------------------------------------------------------
L = 2;
M = 0;
K = 2;
LMK = [L,M,K];
lambda = 2*rand()-1 + 1i*(2*rand()-1);

numtorque = [ -sqrt(3/2)*Z.*(real(lambda)*Y + imag(lambda)*X);
              -sqrt(3/2)*Z.*(real(lambda)*X - imag(lambda)*Y);
               sqrt(3/2)*(2*real(lambda)*X.*Y + imag(lambda)*(X.^2-Y.^2))];

errl20p2val = abs(...
                  runprivate('stochtraj_calcanistorque',LMK,[real(lambda),imag(lambda)],q)-...
                  numtorque...
                 );

errl20p2 = errl20p2val>tol;

% Check for differences between numerical and analytic results
% -------------------------------------------------------------------------

err = any(errl100(:))||any(errl200(:)) ...
      ||any(errl10p1(:))||any(errl20p1(:))||any(errl20p2(:)) ...
      ||any(errl1p10(:))||any(errl1p1p1(:))||any(errl1m1p1(:));

data = [];



end
