% multimatmult  Matlab gateway function to 3D matrix multiplication C 
%               function multimatmult_
%   
%   A = rand(3,3,100); B = rand(3,3,100);
%   C = multimatmult(A,B);
%
%   A = rand(3,3,100) + 1i*rand(3,3,100);
%   B = rand(3,3,100) + 1i*rand(3,3,100);
%   C = multimatmult(A,B);

function C = multimatmult(A,B)

if ismatrix(A)||ismatrix(B)
  mode = '2D';
else
  mode = 'ND';
end

if any(size(A)~=size(B))
  error('Size of input arrays must be equal.')
end

iscomplex = ~isreal(A) || ~isreal(B);

switch mode
  case '2D'
    C = A*B;
  case 'ND'
    switch iscomplex
      case 0
        C = multimatmult_(A,B);
      case 1
        ArBr = multimatmult_(real(A),real(B));
        AiBi = multimatmult_(imag(A),imag(B));
        ArpAiBrpBi = multimatmult_(real(A)+imag(A),real(B)+imag(B));

        C = complex(ArBr - AiBi, ArpAiBrpBi - ArBr - AiBi);
%         C = (ArBr - AiBi) + 1i*(ArpAiBrpBi - ArBr - AiBi);
    end
end

end