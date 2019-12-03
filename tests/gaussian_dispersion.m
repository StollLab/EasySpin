function [err,data] = test(opt,olddata)

% Check Gaussian dispersion lineshape against explicit expression
%=================================================================

% Line parameters
x0 = 2;
fwhm = 0.5;
gamma = fwhm/sqrt(2*log(2)); % distance between inflection points

% Explicitly calculate lineshapes
x = x0 + linspace(-1,1,1e4)*2*gamma;
k = (x-x0)/gamma;

y0 = sqrt(2/pi)/gamma*exp(-2*k.^2);
y1 = sqrt(2/pi)/gamma*2/sqrt(pi)*dawsonF(k*sqrt(2));

% Call gaussian() to get lineshapes
[y2,y3] = gaussian(x,x0,fwhm,0);

if opt.Display
  shift = 0.0;
  plot(x,y0,x,y1,x,y2+shift,x,y3+shift);
  legend('abs explicit','disp explicit','abs gaussian','disp gaussian');
  legend boxoff
end

% Compare
err = ~areequal(y0,y2,1e-10,'rel') || ~areequal(y1,y3,1e-10,'rel');

data =[];

function y = dawsonF(x)
y = real(sqrt(pi)/2i*(faddeeva_(x,38)-exp(-x.^2)));
return

function w = faddeeva_(z,N)

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
