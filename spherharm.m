% spherharm    Spherical and tesseral harmonics
%
%   y = spherharm(L,M,theta,phi)
%   y = spherharm(L,M,theta,phi,'c')
%   y = spherharm(L,M,theta,phi,'s')
%
%   Computes spherical harmonic Y_L,M(theta,phi)
%   with L>=0 and |M|<L. theta and phi can be
%   scalars or arrays (of the same size). If they
%   are arrays, the spherical harmonic is computed
%   element-by-element.
%
%   If the options 'c' or 's' are included, real spherical
%   harmonics called tesseral harmonics (corresponding to
%   the commonly used orbitals in chemistry) are returned.
%   'c' specifies the function containing cos(m*phi), and
%   's' specifies the function containing sin(m*phi).
%   For the tesseral harmonics, M  is restricted
%   to non-negative values. For M=0, the 'c' and the 's'
%   tesseral harmonic and the spherical harmonic are identical.
%
%   The spherical harmonics include the Condon-Shortley
%   phase convention. The tesseral harmonics ('c'
%   and 's') don't include it.

function y = spherharm(L,M,theta,phi,Type);

switch (nargin)
case 0, help(mfilename); return;
case 4,Type = 'e';
case 5, % nop
otherwise, error('Wrong number of parameters!');
end

if (L<0) | (abs(M)>L)
  error('Wrong integer parameters for spherical harmonic.');
end

if strfind(Type,'cs')
  if (M<0)
    error('For tesseral harmonics, M must be non-negative.');
  end
end

% Normalisation constant
p = prod(L-abs(M)+1:L+abs(M));
if (M>0), p = 1/p; end
NormFactor = sqrt((2*L+1)/4/pi*p);

y = NormFactor * plegendre(L,M,cos(theta));

switch (Type)
case 'e', y = (-1)^M  * exp(1i*M*phi) .* y; % Condon-Shortley phase
case 'c', y = sqrt(2) * cos(M*phi) .* y;
case 's', y = sqrt(2) * sin(M*phi) .* y;
otherwise
  error('Wrong string parameter! Must be ''c'' or ''s''.');
end

return
%========================================================

function test

[x,y,z] = sphere(181);
[phi,the] = vec2ang([x(:) y(:) z(:)].');

for N = 0:20
  for M = 0:N

    val = spherharm(N,M,the,phi,'c');
    val = reshape(val,size(x));
    V = val;

    clf; colormap(flipud(jet(256)));
    %colormap gray
    
    subplot(1,2,1); surf(x,y,z,V);
    axis equal, shading interp
    title(sprintf('%d,%d',N,M));
    xlabel('x'); ylabel('y'); zlabel('z');

    subplot(1,2,2);
    h=surf(x.*abs(V),y.*abs(V),z.*abs(V),V);
    axis equal, shading interp
    light, lighting phong
    set(h,'Ambient',0.7);
    xlabel('x'); ylabel('y'); zlabel('z');
    colorbar
    pause;
  end
end



function test2

[x,y,z] = sphere(181);
[phi,the] = vec2ang([x(:) y(:) z(:)].');

for N = 0:20
  for M = -N:N

    val = spherharm(N,M,the,phi);
    val = reshape(val,size(x));
    V = val;

    clf; colormap(flipud(jet(256)));
    %colormap gray
    
    subplot(1,2,1); surf(x,y,z,angle(V));
    axis equal, shading interp
    title(sprintf('%d,%d',N,M));
    xlabel('x'); ylabel('y'); zlabel('z');

    subplot(1,2,2);
    h=surf(x.*abs(V),y.*abs(V),z.*abs(V),V);
    axis equal, shading interp
    light, lighting phong
    set(h,'Ambient',0.7);
    xlabel('x'); ylabel('y'); zlabel('z');
    colorbar
    pause;
  end
end
