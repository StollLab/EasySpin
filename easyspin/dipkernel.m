% dipkernel   Kernel matrix for DEER
%
%  K = dipkernel(t,r)
%  K = dipkernel(t,r,[g1 g2])
%
%  Provides the dipolar kernel matrix for DEER spectroscopy, assuming a
%  complete powder average and no orientation selection.
%
%  Inputs:
%    t ...       array of time values (microseconds)
%    r ...       array of distance values (nanometers)
%    [gA gB] ... the g values of the two spins, assumed 2.002319 if not
%                given
%
%  Output:
%    K ... kernel matrix where K(i,j) is K(t(i),r(j))
%
%  Example: Calculate the DEER signal from a distance distribution P(r)
%
%  r = linspace(0,6,1001);   % distance range
%  P = gaussian(r,3.5,0.3);  % distance distribution
%  t = linspace(-0.2,4,301); % time range
%  K = dipkernel(t,r);       % dipolar kernel
%  V = K*P(:);               % dipolar signal averaged over distribution
%  plot(t,V);
%  

function K = dipkernel(t,r,gAB)

if nargin==0 && nargout==0
  help(mfilename);
end

if nargin<2
  error('At least two inputs are required: t and r.')
end
if nargin>3
  error('No more than three inputs (t,r,g12) are accepted.')
end

if nargin>2
  if numel(gAB)~=2
    error('Third input (g12) must be an array with two elements, [g1 g2].');
  end
  gA = gAB(1);
  gB = gAB(2);
else
  gA = gfree;
  gB = gfree;
end

t = t(:);
r = r(:);

D = (mu0/4/pi)*bmagn^2*gA*gB;
D = D/planck*1e21;  % J m^3 -> MHz nm^3
wdd = 2*pi*D./r.^3;

ph = wdd.'.*abs(t);
kappa = sqrt(6*ph/pi);

C = fresnelC(kappa)./kappa;
S = fresnelS(kappa)./kappa;
K = cos(ph).*C + sin(ph).*S;

% Replace NaN values (at time zero) with 1
K(isnan(K)) = 1;

% Include dr
if numel(r)>1
  dr = r(2)-r(1);
  if abs(dr/mean(diff(r))-1)>1e-8
    error('Distance vector must have equidistant increments.')
  end
  K = dr*K;
end

end
