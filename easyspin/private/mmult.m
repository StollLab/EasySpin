% mmult  
%
%   Matlab gateway function to 3D matrix multiplication C function
%   
%   A = rand(3,3,100); B = rand(3,3,100);
%   C = mmult(A,B,'real');
%
%   A = rand(3,3,100) + 1i*rand(3,3,100);
%   B = rand(3,3,100) + 1i*rand(3,3,100);
%   C = mmult(A,B,'complex');
%

function C = mmult(A,B,complexity)

if ndims(A)<3||ndims(B)<3%size(A,3)<2 || size(B,3)<2
  error('Do not use this function for 2D arrays. Use the "*" operator instead.')
end

if ~ischar(complexity)
  error('Complexity flag must be a character array.')
end

switch complexity
  case 'real'
    C = mmx_simple(A,B);
  case 'complex'
    ArBr = mmx_simple(real(A),real(B));
    AiBi = mmx_simple(imag(A),imag(B));
    ArpAiBrpBi = mmx_simple(real(A)+imag(A),real(B)+imag(B));

    C = (ArBr - AiBi) + 1i*(ArpAiBrpBi - ArBr - AiBi);
  otherwise
    error('Complexity flag must be equal to "real" or "complex".')
end

end