%  autocorrfft  Calculate autocorrelation function using the FFT.
%
%  AutoCorr = autocorrfft(y);
%  AutoCorr = autocorrfft(y, dim);
%  AutoCorr = autocorrfft(y, dim, normalized);
%  AutoCorr = autocorrfft(y, dim, normalized, vector);
%  AutoCorr = autocorrfft(y, dim, normalized, vector, centered);
%
%  Input:
%     y              array of data
%     dim            integer
%                      Dimension over which to perform FFT, first
%                      non-singleton dimension by default
%     normalized     1: normalize by the variance (default)
%                    0: no normalization
%     vector         integer
%                      Dimension over which to sum as components of a
%                      vector (scalar behavior by default)
%     centered       1: subtract mean from data (default)
%                    0: leave data as is
%
%  Output:
%     AutoCorr       autocorrelation of y

function AutoCorr = autocorrfft(varargin)

switch nargin
  case 1
    % normalized output by default
    y = varargin{1};
    dim = [];
    normalized = 1;
    vector = 0;
    centered = 1;
  case 2
    y = varargin{1};
    dim = varargin{2};
    normalized = 1;
    vector = 0;
    centered = 1;
  case 3
    y = varargin{1};
    dim = varargin{2};
    normalized = varargin{3};
    vector = 0;
    centered = 1;
  case 4
    y = varargin{1};
    dim = varargin{2};
    normalized = varargin{3};
    vector = varargin{4};
    centered = 1;
  case 5
    y = varargin{1};
    dim = varargin{2};
    normalized = varargin{3};
    vector = varargin{4};
    centered = varargin{5};
  otherwise
    error('Wrong number of input arguments.')
end

Dimsy = ndims(y);

if isempty(dim)
  % determine first non-singleton dimension
  sizey = size(y);
  idx = find(sizey>1);
  if ~isempty(idx)
    dim = idx(1);
  else
    error('No non-singleton dimensions were detected for input.')
  end
end

% center the data
if centered
  y = bsxfun(@minus, y, mean(y, dim));
end


N = size(y, dim);
F = fft(y, 2*N, dim);
r = ifft(F.*conj(F), [], dim);

% select only the first half of the FFT dimension
idx = cell(1, Dimsy);
idx(:) = {':'};
idx(dim) = {1:N};
AutoCorr = real(r(idx{:}));

if vector
  % sum along dimension given by value of vector argument
  AutoCorr = sum(AutoCorr, vector);
end

n = N*ones(1, N) - [1:N] + 1;
sizen = ones(1,Dimsy);
sizen(dim) = N;
n = reshape(n,sizen);
AutoCorr = AutoCorr./n;

if normalized
  % divide by first value of the FFT dimension
  idx = cell(1, Dimsy);
  idx(:) = {':'};
  idx{dim} = 1;
  AutoCorr = AutoCorr./AutoCorr(idx{:});
end

end