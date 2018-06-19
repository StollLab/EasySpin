% multimatmult  Matlab gateway function to 3D matrix multiplication C function
%   
%   A = rand(3,3,100); B = rand(3,3,100);
%   C = multimatmult(A,B);
%
%   A = rand(3,3,100) + 1i*rand(3,3,100);
%   B = rand(3,3,100) + 1i*rand(3,3,100);
%   C = multimatmult(A,B);
%

function C = multimatmult(A,B)

if ndims(A)<3||ndims(B)<3%size(A,3)<2 || size(B,3)<2
  mode = '2D';
%   error('Do not use this function for 2D arrays. Use the "*" operator instead.')
else
  mode = 'ND';
end

if any(size(A)~=size(B))
  error('Size of input arrays must be equal.')
end

if isreal(A) && isreal(B)
  iscomplex = 0;
else
  iscomplex = 1;
end

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