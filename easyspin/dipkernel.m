% dipkernel   Kernel matrix for DEER
%
%  K = dipkernel(t,r)
%  K = dipkernel(t,r,[g1 g2])
%
%  Provides the kernel matrix for DEER spectroscopy, a complete powder average,
%  and no orientation selection.
%
%  Inputs:
%    t ...       array of time values (microseconds)
%    r ...       array of distance values (nanometers)
%    [g1 g2] ... the g values of the two spins, assumed 2.002319 if not
%                given
%
%  Output:
%    K ... kernel matrix where K(i,j) is K(t(i),r(j))
%

function K = dipkernel(t,r,g12)

if nargin==0, help(mfilename); end

if nargin<2
  error('At least two inputs are required: t and r.')
end
if nargin>3
  error('No more than three inputs (t,r,g12) are accepted.')
end

if nargin>2
  if numel(g12)~=2
    error('Third input (g12) must be an array with two elements, [g1 g2].');
  end
  g1 = g12(1);
  g2 = g12(2);
else
  g1 = gfree;
  g2 = gfree;
end

t = t(:);
r = r(:);

D = (mu0/4/pi)*bmagn^2*g1*g2;
D = D/planck*1e21;  % J m^3 -> MHz nm^3
wdd = 2*pi*D./r.^3;

ph = wdd.'.*abs(t);
kappa = sqrt(6*ph/pi);

C = fresnelC(kappa)./kappa;
S = fresnelS(kappa)./kappa;
K = cos(ph).*C + sin(ph).*S;

% Replace NaN values (at time zero) with 1
K(isnan(K)) = 1;

end
