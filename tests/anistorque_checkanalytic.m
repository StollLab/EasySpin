function [err,data] = test(opt,olddata)
% Check that wignerdquat agrees with analytical D-matrix elements from 
% tables for L = 1,2

rng(1)

tol = 5*eps;
len = 2;

beta  =   pi*(2*rand(1,len)-1);
alpha = 2*pi*(2*rand(1,len)-1);
gamma = 2*pi*(2*rand(1,len)-1);

q = euler2quat(alpha,beta,gamma);

A = q(1,:) - 1i*q(4,:);
B = -q(3,:) - 1i*q(2,:);
Ast = q(1,:) + 1i*q(4,:);
Bst = -q(3,:) + 1i*q(2,:);

X = A.*Bst + Ast.*B;
Y = -1i*(A.*Bst - Ast.*B);
Z = A.*Ast - B.*Bst;

% Note that there are typos in Table 3 of Lynden-Bell, Stone, Mol. Sim. 1989,
% here is the correct matrix, which is actually orthogonal:
%
% R = [ 1/2*(A^2-B^2+Ast^2-Bst^2), -1i/2*(A^2+B^2-Ast^2-Bst^2),      -A*B-Ast*Bst;
%      1i/2*(A^2-B^2-Ast^2+Bst^2),   1/2*(A^2+B^2+Ast^2+Bst^2), -1i*(A*B-Ast*Bst);
%                     A*Bst+Ast*B,           -1i*(A*Bst-Ast*B),      A*Ast-B*Bst];

% Test c10 against Eq. C2
L = 1;
M = 0;
K = 0;
LMK = [L,M,K];
c10 = 1;
lambda = c10 + 1i*c10;
errc10val = abs(anistorque(LMK,[c10,c10],q)-[      -c10*Y;
                                                    c10*X;
                                             zeros(1,len)]);
% errc10val = abs(anistorque(LMK,[c10,c10],q)-[-1i/sqrt(2)*(lambda*(wigD(L,M,K+1,q)+wigD(L,M,K-1,q))+conj(lambda)*(wigD(L,-M,-K+1,q)+wigD(L,-M,-K-1,q)));
%                                               -1/sqrt(2)*(lambda*(wigD(L,M,K+1,q)-wigD(L,M,K-1,q))+conj(lambda)*(wigD(L,-M,-K+1,q)-wigD(L,-M,-K-1,q)));
%                                                zeros(1,len)]/2);
errc10 = errc10val>tol;

% Test c1p1
L = 1;
M = 0;
K = 1;
LMK = [L,M,K];
c1p1 = 1;
lambda = c1p1 + 1i*c1p1;
% errc1p1val = abs(anistorque(LMK,[c1p1,c1p1],q)-[ -1i/sqrt(2)*(lambda*wigD(L,M,K-1,q)...
%                                                               -conj(lambda)*wigD(L,-M,-K+1,q));
%                                                    1/sqrt(2)*(lambda*wigD(L,M,K-1,q)...
%                                                               +conj(lambda)*wigD(L,-M,-K+1,q));
%                                                            -1i*(lambda*wigD(L,M,K,q)...
%                                                                +conj(lambda)*wigD(L,-M,-K,q))]/2);
errc1p1val = abs(anistorque(LMK,[c1p1,c1p1],q)-[                   1/sqrt(2)*imag(lambda)*Z;
                                                                   1/sqrt(2)*real(lambda)*Z;
                                                 -1/sqrt(2)*(real(lambda)*Y+imag(lambda)*X)]);
errc1p1 = errc1p1val>tol;

% anistorque(LMK,[c1p1,c1p1],q)
% 1/sqrt(2)*imag(lambda)*Z

% % % Test c1m1
% L = 1;
% M = 0;
% K = -1;
% LMK = [L,M,K];
% c1m1 = 1;
% lambda = c1m1 + 1i*c1m1;
% % errc1m1val = abs(anistorque(Q,LMK,[c1m1,c1m1])-[      -1i/sqrt(2)*c1m1*Z;
% %                                                           -1/sqrt(2)*c1m1*Z;
% %                                    1i/sqrt(2)*c1m1*(X-1i*Y)]/2;
% % errc1m1val = abs(anistorque(LMK,[c1m1,c1m1],q)-[ -1i/sqrt(2)*(lambda*wigD(L,M,K+1,q)...
% %                                                      -conj(lambda)*wigD(L,-M,-K-1,q));
% %                                                   -1/sqrt(2)*(lambda*wigD(L,M,K+1,q)...
% %                                                      +conj(lambda)*wigD(L,-M,-K-1,q));
% %                                                             1i*(lambda*wigD(L,M,K,q)...
% %                                                        +conj(lambda)*wigD(L,-M,-K,q))]/2);
% errc1m1val = abs(anistorque(LMK,[c1m1,c1m1],q)-[                   1/sqrt(2)*imag(lambda)*Z;
%                                                                    -1/sqrt(2)*real(lambda)*Z;
%                                                  -1/sqrt(2)*(real(lambda)*Y+imag(lambda)*X)])
% errc1m1 = errc1m1val>tol;


% Test c20 against Eq. C8
L = 2;
M = 0;
K = 0;
LMK = [L,M,K];
c20 = 1;
lambda = c20 - 1i*c20;
% errc20val = abs(anistorque(Q,LMK,[c20,c20])-[   -3*c20*Y.*Z;
%                                                  3*c20*X.*Z;
%                                               zeros(1,len)])>tol;
errc20val = abs(anistorque(LMK,[c20,c20],q)-[-1i*sqrt(3/2)*(lambda*(wigD(L,M,K+1,q)+wigD(L,M,K-1,q))...
                                                   +conj(lambda)*(wigD(L,-M,-K+1,q)+wigD(L,-M,-K-1,q)));
                                                -sqrt(3/2)*(lambda*(wigD(L,M,K+1,q)-wigD(L,M,K-1,q))...
                                                   +conj(lambda)*(wigD(L,-M,-K+1,q)-wigD(L,-M,-K-1,q)));
                                                                              zeros(1,len)]/2);
errc20 = errc20val>tol;

% Test c2p1
L = 2;
M = 0;
K = 1;
LMK = [L,M,K];
c2p1 = 1;
lambda = c2p1 + 1i*c2p1;
% errc2p1 = abs(anistorque(Q,LMK,[c2p1,c2p1])-[-1i*sqrt(6)/4*c2p1*((X+1i*Y).^2+(3*Z.^2-1));
%                                     -1*sqrt(6)/4*c2p1*((X+1i*Y).^2-(3*Z.^2-1));
%                                                 1i*sqrt(3/2)*c2p1*(X+1i*Y).*Z])>tol;
% errc2p1val = abs(anistorque(LMK,[c2p1,c2p1],q)-[-1i/2*(lambda*(2*wigD(L,M,K+1,q)+sqrt(6)*wigD(L,M,K-1,q))...
%                                                    -conj(lambda)*(sqrt(6)*wigD(L,-M,-K+1,q)+2*wigD(L,-M,-K-1,q)));
%                                                 -1/2*(lambda*(2*wigD(L,M,K+1,q)-sqrt(6)*wigD(L,M,K-1,q))...
%                                                    -conj(lambda)*(sqrt(6)*wigD(L,-M,-K+1,q)-2*wigD(L,-M,-K-1,q)));
%                                                             -1i*(lambda*wigD(L,M,K,q)+conj(lambda)*wigD(L,-M,-K,q))]/2);
errc2p1val = abs(anistorque(LMK,[c2p1,c2p1],q)-[ sqrt(3/2)/2*(2*real(lambda)*X.*Y+imag(lambda)*(3*Z.^2+X.^2-Y.^2-1));
                                                 sqrt(3/2)/2*(2*imag(lambda)*X.*Y+real(lambda)*(3*Z.^2-X.^2+Y.^2-1));
                                                 -sqrt(3/2)*(imag(lambda)*X+real(lambda)*Y).*Z]);

errc2p1 = errc2p1val>tol;

% Test c2p2 against Eq. C10
L = 2;
M = 0;
K = 2;
LMK = [L,M,K];
c2p2 = 1;
lambda = c2p2 + 1i*c2p2;

% errc2p2val = abs(anistorque(LMK,[c1p1,c1p1],q)-[ -1i*(lambda*wigD(L,M,K-1,q)...
%                                                   +conj(lambda)*wigD(L,-M,-K+1,q));
%                                                  (lambda*wigD(L,M,K-1,q)...
%                                                   -conj(lambda)*wigD(L,-M,-K+1,q));
%                                                  -2i*(lambda*wigD(L,M,K,q)...
%                                                   -conj(lambda)*wigD(L,-M,-K,q))]/2);

errc2p2val = abs(anistorque(LMK,[c2p2,c2p2],q)-[           -sqrt(3/2)*Z.*(real(lambda)*Y + imag(lambda)*X);
                                                           -sqrt(3/2)*Z.*(real(lambda)*X - imag(lambda)*Y);
                                                sqrt(3/2)*(2*real(lambda)*X.*Y + imag(lambda)*(X.^2-Y.^2))]);
% anistorque(LMK,[c2p2,c2p2],q)
% -sqrt(6)*c2p2*Y.*Z
% -sqrt(6)*c2p2*X.*Z


errc2p2 = errc2p2val>tol;

err = any(errc10(:))||any(errc1p1(:)) ...
      ||any(errc20(:))||any(errc2p1(:))||any(errc2p2(:));

data = [];


%% Helper functions
function Del = wigD(L,M,K,q)
% Helper function for calling L=1 Wigner D-matrix elements

A = q(1,:) - 1i*q(4,:);
B = -q(3,:) - 1i*q(2,:);
Ast = q(1,:) + 1i*q(4,:);
Bst = -q(3,:) + 1i*q(2,:);

X = A.*Bst + Ast.*B;
Y = -1i*(A.*Bst - Ast.*B);
Z = A.*Ast - B.*Bst;

switch L
  case 1
    D = zeros(3,3,size(A,2));
    
    D(1,1,:) = A.^2;
    D(1,2,:) = sqrt(2)*A.*B;
    D(1,3,:) = B.^2;
    D(2,1,:) = -sqrt(2)*A.*Bst;
    D(2,2,:) = Z;
    D(2,3,:) = sqrt(2)*Ast.*B;
    D(3,1,:) = Bst.^2;
    D(3,2,:) = -sqrt(2).*Ast.*Bst;
    D(3,3,:) = Ast.^2;

    Del = reshape(D(2-M,2-K,:),[1,size(q,2)]);
  case 2
    D = zeros(5,5,size(A,2));
    
    D(1,1,:) = A.^4;
    D(1,2,:) = 2*A.^3.*B;
    D(1,3,:) = sqrt(6)*A.^2.*B.^2;
    D(1,4,:) = 2*A.*B.^3;
    D(1,5,:) = B.^4;
    D(2,1,:) = -2*A.^3.*Bst;
    D(2,2,:) = A.^2.*(2*Z-1);
    D(2,3,:) = sqrt(6)*A.*B.*Z;
    D(2,4,:) = B.^2.*(2*Z+1);
    D(2,5,:) = 2*Ast.*B.^3;
    D(3,1,:) = sqrt(6)*A.^2.*Bst.^2;
    D(3,2,:) = -sqrt(6)*A.*Bst.*Z;
    D(3,3,:) = 1/2*(3*Z.^2-1);
    D(3,4,:) = sqrt(6)*Ast.*B.*Z;
    D(3,5,:) = sqrt(6)*Ast.^2.*B.^2;
    D(4,1,:) = -2*A.*Bst.^3;
    D(4,2,:) = Bst.^2.*(2*Z+1);
    D(4,3,:) = -sqrt(6)*Ast.*Bst.*Z;
    D(4,4,:) = Ast.^2.*(2*Z-1);
    D(4,5,:) = 2*Ast.^3.*B;
    D(5,1,:) = Bst.^4;
    D(5,2,:) = -2*Ast.*Bst.^3;
    D(5,3,:) = sqrt(6)*Ast.^2.*Bst.^2;
    D(5,4,:) = -2*Ast.^3.*Bst;
    D(5,5,:) = Ast.^4;
    
    Del = reshape(D(3-M,3-K,:),[1,size(q,2)]);
end

end

function C = Cpm(L, K, pm)
% Eq. 60 in [1]

if     pm=='+', C = sqrt(L*(L+1)-K*(K+1));
elseif pm=='-', C = sqrt(L*(L+1)-K*(K-1)); end

end

end
