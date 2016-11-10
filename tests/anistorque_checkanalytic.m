function [err,data] = test(opt,olddata)
% Check that wignerdquat agrees with analytical D-matrix elements from 
% tables for j = 1,2

tol = 5*eps;
len = 10;

beta = pi*(2*rand(1,len)-1);
alpha = 2*pi*(2*rand(1,len)-1);
gamma = 2*pi*(2*rand(1,len)-1);

q = euler2quat(alpha,beta,gamma);
Q(1,1,:) =  q(1,:) - 1i*q(4,:);
Q(1,2,:) = -q(3,:) - 1i*q(2,:);
Q(2,1,:) =  q(3,:) - 1i*q(2,:);
Q(2,2,:) =  q(1,:) + 1i*q(4,:);

A = q(1,:) - 1i*q(4,:);
B = -q(3,:) - 1i*q(2,:);
Ast = q(1,:) + 1i*q(4,:);
Bst = -q(3,:) + 1i*q(2,:);

% Note that there are typos in Table 3 of Lynden-Bell, Stone, Mol. Sim. 1989,
% here is the correct matrix, which is actually orthogonal:
%
% R = [ 1/2*(A^2-B^2+Ast^2-Bst^2), -1i/2*(A^2+B^2-Ast^2-Bst^2),      -A*B-Ast*Bst;
%      1i/2*(A^2-B^2-Ast^2+Bst^2),   1/2*(A^2+B^2+Ast^2+Bst^2), -1i*(A*B-Ast*Bst);
%                     A*Bst+Ast*B,           -1i*(A*Bst-Ast*B),      A*Ast-B*Bst];

X = A.*Bst + Ast.*B;
Y = -1i*(A.*Bst - Ast.*B);
Z = A.*Ast - B.*Bst;

% Test c10 against Eq. C2
c10 = 2*rand()-1;
coefs.c10 = c10;
errc10 = abs(anistorque(Q,coefs)-[-c10*Y;
                                   c10*X;
                                  zeros(size(Z))])>tol;
coefs = rmfield(coefs,'c10');

% Test c1p1
c1p1 = 2*rand()-1;
coefs.c1p1 = c1p1;
% errc1p1 = abs(anistorque(Q,coefs)-[-1i/sqrt(2)*c1p1*Z;
%                                      1/sqrt(2)*c1p1*Z;
%                                    1i/sqrt(2)*c1p1*(X+1i*Y)])>tol;
errval = abs(anistorque(Q,coefs)-[-1i/sqrt(2)*c1p1*Z;
                                     1/sqrt(2)*c1p1*Z;
                                   1i/sqrt(2)*c1p1*(X+1i*Y)]);
errc1p1 = errval>tol;
coefs = rmfield(coefs,'c1p1');

% % Test c1m1
c1m1 = 2*rand()-1;
coefs.c1m1 = c1m1;
errc1m1 = abs(anistorque(Q,coefs)-[-1i/sqrt(2)*c1m1*Z;
                                    -1/sqrt(2)*c1m1*Z;
                                   1i/sqrt(2)*c1m1*(X-1i*Y)])>tol;
coefs = rmfield(coefs,'c1m1');

% Test c20 against Eq. C8
c20 = 2*rand()-1;
coefs.c20 = c20;
errc20 = abs(anistorque(Q,coefs)-[-3*c20*Y.*Z;
                                   3*c20*X.*Z;
                                  zeros(size(Z))])>tol;
coefs = rmfield(coefs,'c20');

% Test c2p1
c2p1 = 2*rand()-1;
coefs.c2p1 = c2p1;
errc2p1 = abs(anistorque(Q,coefs)-[-1i*sqrt(6)/4*c2p1*((X+1i*Y).^2+(3*Z.^2-1));
                                    -1*sqrt(6)/4*c2p1*((X+1i*Y).^2-(3*Z.^2-1));
                                               1i*sqrt(3/2)*c2p1*(X+1i*Y).*Z])>tol;
coefs = rmfield(coefs,'c2p1');

% Test c2p2 against Eq. C10
c2p2 = 2*rand()-1;
coefs.c2p2 = c2p2;
errc2p2 = abs(anistorque(Q,coefs)-[1i*sqrt(6)*c2p2*A.*Bst.*Z;
                                   -1*sqrt(6)*c2p2*A.*Bst.*Z;
                                 -2i*sqrt(6)*c2p2*A.^2.*Bst.^2])>tol;                      

err = any(errc10(:))||any(errc1p1(:))||any(errc1m1(:)) ...
      ||any(errc20(:))||any(errc2p1(:))||any(errc2p2(:));

data = [];
