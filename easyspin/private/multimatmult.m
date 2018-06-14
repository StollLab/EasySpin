% multimatmult  Matlab gateway function to 3D matrix multiplication C function
%   
%   A = rand(3,3,100); B = rand(3,3,100);
%   C = multimatmult(A,B,'real');
%
%   A = rand(3,3,100) + 1i*rand(3,3,100);
%   B = rand(3,3,100) + 1i*rand(3,3,100);
%   C = multimatmult(A,B,'complex');
%

function C = multimatmult(A,B,complexity)

if ndims(A)<3||ndims(B)<3%size(A,3)<2 || size(B,3)<2
  mode = '2D';
%   error('Do not use this function for 2D arrays. Use the "*" operator instead.')
else
  mode = 'ND';
end

if any(size(A)~=size(B))
  error('Size of input arrays must be equal.')
end

if exist('complexity', 'var')==0
  complexity = 'complex';
end

if ~ischar(complexity)
  error('Complexity flag must be a character array.')
end

switch mode
  case '2D'
    C = A*B;
  case 'ND'
    switch complexity
      case 'real'
        C = multimatmult_(A,B);
      case 'complex'
        ArBr = multimatmult_(real(A),real(B));
        AiBi = multimatmult_(imag(A),imag(B));
        ArpAiBrpBi = multimatmult_(real(A)+imag(A),real(B)+imag(B));

        C = (ArBr - AiBi) + 1i*(ArpAiBrpBi - ArBr - AiBi);
      otherwise
        error('Complexity flag must be equal to "real" or "complex".')
    end
end

end