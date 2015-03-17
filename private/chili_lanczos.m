function [alpha,beta,specerr] = chili_lanczos(A,Vector,z,Options)

nDim = length(Vector);
alpha = zeros(1,nDim);
beta = zeros(1,nDim);
specerr = NaN(1,nDim);

oldspec = inf;
zTest = linspace(z(1),z(end),201);
RecomputationInterval = min(10,ceil(nDim/20));

% Lanczos initializations
x = Vector;
x = x/sqrt(x.'*x);
y = 0;

LentzMethod = Options.Lentz;

if LentzMethod
  %======================================
  tiny = 1e-30;
  spec = tiny;
  C = spec;
  D = 0;
  %======================================
end

% Lanczos iterations
for k = 1:nDim
  if LentzMethod
    %======================================
    Cold = C;
    Dold = D;
    specold = spec;
    %======================================
  end

  % Lanczos
  y = A*x + y;
  alpha(k) = x.'*y;
  y = y - alpha(k)*x;
  beta(k) = sqrt(y.'*y);
  tempc = x;
  x = y/beta(k);
  y = -beta(k)*tempc;

  % Continued fraction: next convergent
  if LentzMethod
    %=====================================
    a = z + alpha(k);
    if (k==1), b = 1; else b = -beta(k-1)^2; end
    C = a + b./Cold;
    D = 1./(a + b.*Dold);
    delta = C.*D;
    spec = specold.*delta;
    %specerr(k) = max(abs(delta-1));
    if ~mod(k,RecomputationInterval)
      specerr(k) = norm(delta-1,inf);
      if (specerr(k)<Options.Threshold), break; end
    end
    %=====================================
  else
    if ~mod(k,RecomputationInterval)
      spec = chili_contfracspec(zTest,alpha,beta,k);
      respec = real(spec);
      specerr(k) = max(abs(respec-oldspec))/max(respec);
      if (specerr(k)<Options.Threshold), break; end
      oldspec = respec;
    end
  end

  if (beta(k)==0)
    specerr(k) = 0;
    break
  end
  
end

alpha = alpha(1:k);
beta = beta(1:k);
specerr = specerr(1:k);

return
