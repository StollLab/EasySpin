% faddeeva  Faddeeva function
%
%   y = faddeeva(z)
%   y = faddeeva(z,N)
%
%   Returns the Faddeeva function w(z) evaluated at the points in z.
%   N determines the number of terms in the expansion.

function w = faddeeva(z,N)

% The Faddeeva function is
%
%      w(z) = exp(-x^2)erfc(-iz)
%           = exp(-z^2)(1+2i/sqrt(pi)*\int_0^z e^(t^2) dt)
% 
% Below an implementation of Weideman's algorithm. See
% Weideman, J. A. C., Computation of the complex error function.
% SIAM J. Numer. Anal. 31 (1994), no. 5, 1497-1518.

if (nargin<2), N = 38; end

M = 2*N;
% M2 = no. of sampling points
M2 = 2*M;                                    
k = (-M+1:1:M-1).';

% Optimal choice of L
L = sqrt(N/sqrt(2));

% Define variables theta and t
theta = k*pi/M;
t = L*tan(theta/2);

% Function to be transformed
f = exp(-t.^2).*(L^2+t.^2);
f = [0; f];

a = real(fft(fftshift(f)))/M2;               % Coefficients of transform
a = flipud(a(2:N+1));                        % Reorder coefficients

Z = (L+1i*z)./(L-1i*z);
p = polyval(a,Z);    % Polynomial evaluation

w = 2*p./(L-1i*z).^2+(1/sqrt(pi))./(L-1i*z); % Evaluate w(z)
