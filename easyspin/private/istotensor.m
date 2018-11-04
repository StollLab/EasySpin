% istotensor  Irreducbile spherical tensor operator from 2 vector operators
%
%   [T0,T1,T2] = istotensor(a,b)
%
% istotensor converts two cartesian vector operators (or, more specifically,
% their outer product b*a') into three irreducible spherical tensors
% T^(l)(a,b) of rank l=0, 1, and 2.
%
% Examples of vectors/vector operators are the magnetic field, electron spin
% vector operators, or nuclear spin vector operators.
%
% Input:
%    a  ... 3-element array or cell array for the first vector/vector operator
%    b  ... 3-element array or cell array for the second vector/vector operator
%
% Output:
%    T0  ... rank-0 irreducbile spherical tensor (scalar)
%    T1  ... rank-1 irreducible spherical tensor (3x1 cell array)
%    T2  ... rank-2 irreducible spherical tensor (5x1 cell array)
%
% Reference:
%   Michael Mehring
%   Principles of High Resolution NMR in Solids, 2nd edition
%   Wiley, 1983
%   Appendix A, starting on p.288

function [T0,T1,T2] = istotensor(a,b)

if iscell(a)
  ax = a{1}; ay = a{2}; az = a{3};
else
  ax = a(1); ay = a(2); az = a(3);
end
if iscell(b)
  bx = b{1}; by = b{2}; bz = b{3};
else
  bx = b(1); by = b(2); bz = b(3);
end
  
T0 = -(1/sqrt(3)) * (ax*bx+ay*by+az*bz);             % (0,0)

T1 = cell(3,1);
T1{1} = -0.5*((ax*bz - az*bx) + 1i*(ay*bz - az*by)); % (1,+1)
T1{2} = -(1i/sqrt(2))*(ay*bx - ax*by);               % (1,0)
T1{3} = -0.5*((ax*bz - az*bx) - 1i*(ay*bz - az*by)); % (1,-1)

T2 = cell(5,1);
T2{1} =  0.5*((ax*bx - ay*by) + 1i*(ax*by + ay*bx)); % (2,+2)
T2{2} = -0.5*((ax*bz + az*bx) + 1i*(ay*bz + az*by)); % (2,+1)
T2{3} =  sqrt(2/3) *  (az*bz - 0.5*(ax*bx + ay*by)); % (2,0)
T2{4} =  0.5*((ax*bz + az*bx) - 1i*(ay*bz + az*by)); % (2,-1)
T2{5} =  0.5*((ax*bx - ay*by) - 1i*(ax*by + ay*bx)); % (2,-2)
  
return
