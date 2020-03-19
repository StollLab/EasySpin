% tensor_sph2cart  Rank-2 cartesian tensor from irreducible spherical tensors
%
%    Tc = tensor_sph2cart(T0,T1,T2)
%
% tensor_sph2cart converts the three irreducible spherical tensors (rank 0, 1,
% and 2) to the corresponding rank-2 cartesian tensor.
%
% Input:
%    T0 ... rank-0 spherical tensor (1 element, real-valued)
%    T1 ... rank-1 spherical tensor (3 elements)
%    T2 ... rank-2 spherical tensor (5 elements)
%
% Output:
%    Tc ... rank-2 cartesian tensor (3x3 array, real-valued)

% Reference:
%   Michael Mehring
%   Principles of High Resolution NMR in Solids, 2nd edition
%   Wiley, 1983
%   Appendix A, starting on p.288

function T = tensor_sph2cart(T0,T1,T2)

% Check inputs
%-------------------------------------------------------------------------------
if nargin~=3
  error('Three inputs (T0, T1, T2) are required.');
end

if ~isnumeric(T0) || numel(T0)~=1
  error('First input (T0) must be a single number.');
end
if ~isreal(T0)
  error('First input (T0) must be a real-valued number.');
end

if ~isnumeric(T1) || numel(T1)~=3
  error('Second input (T1) must be an array with 3 numbers.');
end

if ~isnumeric(T2) || numel(T2)~=5
  error('Third input (T2) must be an array with 5 numbers.');
end

% Calculate cartesian tensor elements
%-------------------------------------------------------------------------------

% Cartesian indices (for better readability)
x = 1;
y = 2;
z = 3;

% Define functions that take q as input (for better readability)
T1 = @(q)T1(2-q);
T2 = @(q)T2(3-q);

T(x,z) =  0.5 *T1(1) + 0.5 *T1(-1) - 0.5 *T2(1) + 0.5 *T2(-1);
T(z,x) = -0.5 *T1(1) - 0.5 *T1(-1) - 0.5 *T2(1) + 0.5 *T2(-1);
T(y,z) = -0.5i*T1(1) + 0.5i*T1(-1) + 0.5i*T2(1) + 0.5i*T2(-1);
T(z,y) =  0.5i*T1(1) - 0.5i*T1(-1) + 0.5i*T2(1) + 0.5i*T2(-1);

T(x,x) = -1/sqrt(3)*T0 + 0.5*T2(2) - 1/sqrt(6)*T2(0) + 0.5*T2(-2);
T(y,y) = -1/sqrt(3)*T0 - 0.5*T2(2) - 1/sqrt(6)*T2(0) - 0.5*T2(-2);
T(z,z) = -1/sqrt(3)*T0 + sqrt(2/3)*T2(0);
T(x,y) =  1i/sqrt(2)*T1(0) - 0.5i*T2(2) + 0.5i*T2(-2);
T(y,x) = -1i/sqrt(2)*T1(0) - 0.5i*T2(2) + 0.5i*T2(-2);

return
