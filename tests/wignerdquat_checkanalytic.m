function [err,data] = test(opt,olddata)
% Check that wignerdquat agrees with analytical results from tables for
% L = 1,2

len = 5;

beta = pi*(2*rand(1,len)-1);
alpha = 2*pi*(2*rand(1,len)-1);
gamma = 2*pi*(2*rand(1,len)-1);

q = euler2quat(alpha,beta,gamma);

A = q(1,:) - 1i*q(4,:);
B = -q(3,:) - 1i*q(2,:);
Ast = q(1,:) + 1i*q(4,:);
Bst = -q(3,:) + 1i*q(2,:);

Z = A.*Ast - B.*Bst;

% D matrix for L=1, ordered such that rows (1,2,3)<->M=(+1,0,-1), and
% cols (1,2,3)<->K=(+1,0,-1)
% D1 = [           A.^2,      sqrt(2)*A.*B,            B.^2;
%       -sqrt(2)*A.*Bst,                Z,   sqrt(2)*Ast.*B;
%                Bst.^2, -sqrt(2).*Ast.*Bst,         Ast.^2];
D1(1,1,:) = A.^2;
D1(1,2,:) = sqrt(2)*A.*B;
D1(1,3,:) = B.^2;
D1(2,1,:) = -sqrt(2)*A.*Bst;
D1(2,2,:) = Z;
D1(2,3,:) = sqrt(2)*Ast.*B;
D1(3,1,:) = Bst.^2;
D1(3,2,:) = -sqrt(2).*Ast.*Bst;
D1(3,3,:) = Ast.^2;

% D matrix for L=2, ordered such that rows (1,2,3,4,5)<->M=(2,1,0,-1,-2),
% and cols (1,2,3)<->K=(2,1,0)
% D2 = [                A.^4,          2*A.^3.*B,     sqrt(6)*A.^2.*B.^2;
%               -2*A.^3.*Bst,      A.^2.*(2*Z-1),        sqrt(6)*A.*B.*Z;
%       sqrt(6)*A.^2.*Bst.^2, -sqrt(6)*A.*Bst.*Z,         1/2*(3*Z.^2-1);
%               -2*A.*Bst.^3,    Bst.^2.*(2*Z+1),   -sqrt(6)*Ast.*Bst.*Z;
%                     Bst.^4,     -2*Ast.*Bst.^3, sqrt(6)*Ast.^2.*Bst.^2];
% D2(1,1,:) = A.^4;
% D2(1,2,:) = 2*A.^3.*B;
% D2(1,3,:) = sqrt(6)*A.^2.*B.^2;
% D2(2,1,:) = -2*A.^3.*Bst;
% D2(2,2,:) = A.^2.*(2*Z-1);
% D2(2,3,:) = sqrt(6)*A.*B.*Z;
% D2(3,1,:) = sqrt(6)*A.^2.*Bst.^2;
% D2(3,2,:) = -sqrt(6)*A.*Bst.*Z;
% D2(3,3,:) = 1/2*(3*Z.^2-1);
% D2(4,1,:) = -2*A.*Bst.^3;
% D2(4,2,:) = Bst.^2.*(2*Z+1);
% D2(4,3,:) = -sqrt(6)*Ast.*Bst.*Z;
% D2(5,1,:) = Bst.^4;
% D2(5,2,:) = -2*Ast.*Bst.^3;
% D2(5,3,:) = sqrt(6)*Ast.^2.*Bst.^2;
D2(1,1,:) = A.^4;
D2(1,2,:) = 2*A.^3.*B;
D2(1,3,:) = sqrt(6)*A.^2.*B.^2;
D2(1,4,:) = 2*A.*B.^3;
D2(1,5,:) = B.^4;
D2(2,1,:) = -2*A.^3.*Bst;
D2(2,2,:) = A.^2.*(2*Z-1);
D2(2,3,:) = sqrt(6)*A.*B.*Z;
D2(2,4,:) = B.^2.*(2*Z+1);
D2(2,5,:) = 2*Ast.*B.^3;
D2(3,1,:) = sqrt(6)*A.^2.*Bst.^2;
D2(3,2,:) = -sqrt(6)*A.*Bst.*Z;
D2(3,3,:) = 1/2*(3*Z.^2-1);
D2(3,4,:) = sqrt(6)*Ast.*B.*Z;
D2(3,5,:) = sqrt(6)*Ast.^2.*B.^2;
D2(4,1,:) = -2*A.*Bst.^3;
D2(4,2,:) = Bst.^2.*(2*Z+1);
D2(4,3,:) = -sqrt(6)*Ast.*Bst.*Z;
D2(4,4,:) = Ast.^2.*(2*Z-1);
D2(4,5,:) = 2*Ast.^3.*B;
D2(5,1,:) = Bst.^4;
D2(5,2,:) = -2*Ast.*Bst.^3;
D2(5,3,:) = sqrt(6)*Ast.^2.*Bst.^2;
D2(5,4,:) = -2*Ast.^3.*Bst;
D2(5,5,:) = Ast.^4;
                
              
L = 1;
err1 = zeros(3,3,len);

for M=-L:L
  for K=-L:L
    err1(2-M,2-K,:) = abs(wignerdquat(L,M,K,q)-reshape(D1(2-M,2-K,:),[1,len]));
  end
end

L = 2;
err2 = zeros(5,5,len);

for M=-L:L
  for K=-L:L
    err2(3-M,3-K,:) = abs(wignerdquat(L,M,K,q)-reshape(D2(3-M,3-K,:),[1,len]));
  end
end
err = any(err1(:)>1e-15)||any(err2(:)>1e-15);

data = [];
