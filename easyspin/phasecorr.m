% phasecorr  Phase-correct data to minimize imaginary part
%
%  Vph = phasecorr(V)
%  Vph = phasecorr(V,ignoreOffset)
%  [Vph,ph] = phasecorr(___)
%
%  V  ... one- or two-dimensional array of data
%  ignoreOffset... whether to ignore the offset of the imaginary
%      part (default false)
%  ph ... phase correction angle, in radians
%  Vph ... phase-corrected signal, equal to V.*exp(1i*ph)
%
%  This function determines the angle ph such that it minimizes
%  the second moment (if ignoreOffset is set to false) or the
%  second central moment (if ignoreOffset is set to true) of the
%  imaginary part.
%
%  Example:
%  f = 2.2134;  % frequency, MHz
%  t = linspace(0,5,1001);  % time, µs
%  V = cos(2*pi*f*t).*exp(-1i*deg2rad(18));  % signal
%  Vph = phasecorr(V);
%  subplot(2,1,1); plot(t,real(V),t,imag(V));
%  subplot(2,1,2); plot(t,real(Vph),t,imag(Vph));

function varargout = phasecorr(V,ignoreOffset)

% The following determines the phase that minimizes the objective
% function = sum of squares of imaginary part of V*exp(1i*phi)
% This objective function has the analytical form
%
%   (a+b) + (b-a)*cos(2*phi) + c*sin(2*phi)
%    = offset + amp*cos(2*phi-phi0)
%
% where
%    a = sum_k real(V_k)^2 / 2
%    b = sum_k imag(V_k)^2 / 2
%    c = sum_k real(V_k)*imag(V_k)
%
%    offset = a+b
%    amp = sqrt((b-a)^2+c^2)
%    phi0 = atan2(c,b-a)
%
% The objective function has two minima:
%    phi = phi0/2 + pi/2
%    phi = phi0/2 + 3*pi/2

if nargin==0
  help(mfilename);
  return
end

if nargin<2
  ignoreOffset = false;
end

% Determine phase that minimizes imaginary part (with or without offset)
if ignoreOffset
  phimin = fminbnd(@(ph)objfun(V,ph),0,pi);
else
  Vr = real(V);
  Vi = imag(V);
  a = sum(Vr.^2,'all')/2;
  b = sum(Vi.^2,'all')/2;
  c = sum(Vr.*Vi,'all');
  phi0 = atan2(c,b-a);
  phimin = phi0/2 + pi/2;
end

% Determine phase that results in smaller phase shift
phaseshift = @(a) mod(a+pi,2*pi)-pi;
if abs(phaseshift(phimin+pi))<abs(phaseshift(phimin))
  phimin = phimin + pi;
end

% Wrap phase shift angle to (-pi,pi)
phimin = mod(phimin+pi,2*pi)-pi;

V_phased = V*exp(1i*phimin);

switch nargout
  case 0
    if isvector(V)
      plotphasecorr(V,V_phased,phimin);
    else
      error('Plotting of multidimensional data is not implemented.');
    end
    varargout = {};
  case 1
    varargout = {V_phased};
  case 2
    varargout  = {V_phased,phimin};
end

end

% Objective function for phasing
%-----------------------------------------------------------------------------
function f = objfun(V,phi)
Vi = imag(V*exp(1i*phi));
Vic = Vi - mean(Vi,'all');
f = sum(Vic.^2,'all');
end

% Plotting
%-----------------------------------------------------------------------------
function plotphasecorr(V,V_phased,phimin)

n = numel(V);

tiledlayout(2,3,Padding="compact");
nexttile(1,[1,2])
plot(1:n,real(V),1:n,imag(V));
title('original data V')
legend('Re','Im');
legend boxoff

nexttile(3)
polarplot(angle(V),abs(V));

nexttile(4,[1 2])
plot(1:n,real(V_phased),1:n,imag(V_phased))
title(sprintf('phased data V*exp(1i*phase), phase = %0.4f deg',rad2deg(phimin)))
legend('Re','Im');
legend boxoff

nexttile(6)
polarplot(angle(V_phased),abs(V_phased));

end
