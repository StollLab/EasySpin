% randi() generate random integer
%
%  X = randi(imax)
%  X = randi(imax,N,M)
%
%

function X = randi_(varargin)

if verLessThan('matlab','7.12')
  switch nargin
    case 1
      imax = varargin{1};
      N = 1;
      M = 1;
    case 3
      imax = varargin{1};
      N = varargin{2};
      M = varargin{3};
    otherwise
      error('Need one (imax) or three (imax, N, M) input arguments.')
  end

  X = rand(N,M);
  X = ceil(X*imax);
else
  X = randi(varargin{:});
end
