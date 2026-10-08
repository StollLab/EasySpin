% strains_modes  Mode matrix from strain covariance input
%
%   Q = strains_modes(Sys,n)
%
%   Validates Sys.StrainFWHM/Sys.StrainCorr or Sys.StrainModes for n strain
%   parameters and returns the n x m mode matrix Q, with covariance Q*Q.'
%   (in FWHM units). Columns that are zero are dropped. Errors are thrown
%   directly.

function Q = strains_modes(Sys,n)

hasFWHM = isfield(Sys,'StrainFWHM');
hasCorr = isfield(Sys,'StrainCorr');
hasModes = isfield(Sys,'StrainModes');

fields = {'StrainFWHM','StrainCorr','StrainModes'};
for f = 1:numel(fields)
  if isfield(Sys,fields{f}) && isempty(Sys.(fields{f}))
    error('Sys.%s is empty.',fields{f});
  end
end

if hasFWHM==hasModes
  error('With Sys.StrainPars, give either Sys.StrainFWHM (optionally with Sys.StrainCorr) or Sys.StrainModes.');
end
if hasCorr && ~hasFWHM
  error('Sys.StrainCorr can only be used with Sys.StrainFWHM.');
end

if hasFWHM

  fwhm = Sys.StrainFWHM;
  if ~isnumeric(fwhm) || ~isreal(fwhm) || ~isvector(fwhm) || numel(fwhm)~=n || ...
      any(~isfinite(fwhm)) || any(fwhm<0)
    error('Sys.StrainFWHM must contain %d nonnegative real numbers, one per entry in Sys.StrainPars.',n);
  end

  % Correlation matrix
  if hasCorr
    c = Sys.StrainCorr;
    if ~isnumeric(c) || ~isreal(c) || any(~isfinite(c(:)))
      error('Sys.StrainCorr must contain real numbers.');
    end
    nUpper = n*(n-1)/2;
    if n>1 && isvector(c) && numel(c)==nUpper
      % upper triangle in row order (1-2, 1-3, ..., 2-3, ...)
      C = zeros(n);
      C(tril(true(n),-1)) = c;
      C = C + C.' + eye(n);
    elseif isequal(size(c),[n n])
      C = c;
    else
      error('Sys.StrainCorr must be a %dx%d matrix, or contain the %d elements of its upper triangle.',n,n,nUpper);
    end
    if norm(C-C.',1)>1e-12
      error('Sys.StrainCorr must be symmetric.');
    end
    if any(abs(diag(C)-1)>1e-12)
      error('Sys.StrainCorr must have ones on the diagonal.');
    end
    if any(abs(C(:))>1+1e-12)
      error('Sys.StrainCorr elements must be between -1 and +1.');
    end
  else
    C = eye(n);
  end

  % Mode matrix Q = S*V*sqrt(Lambda), with C = V*Lambda*V.'
  [V,L] = eig((C+C.')/2);
  L = diag(L);
  if min(L)<-1e-10
    error('Sys.StrainCorr is not a valid correlation matrix (it is not positive semidefinite).');
  end
  L(L<0) = 0;
  Q = diag(fwhm)*V*diag(sqrt(L));

else

  modes = Sys.StrainModes;
  if ~isnumeric(modes) || ~isreal(modes) || ndims(modes)~=2 || size(modes,2)~=n || any(~isfinite(modes(:)))
    error('Sys.StrainModes must be a real array with %d columns, one per entry in Sys.StrainPars.',n);
  end
  Q = modes.';

end

% Drop zero modes
colNorm = sqrt(sum(Q.^2,1));
Q = Q(:,colNorm>1e-12*max([colNorm 0]));

end
