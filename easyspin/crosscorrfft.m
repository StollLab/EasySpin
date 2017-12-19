%  crosscorrfft  Calculate cross-correlation function using the FFT.
%
%  CrossCorr = crosscorrfft(x, y);
%  CrossCorr = crosscorrfft(x, y, dim);
%  CrossCorr = crosscorrfft(x, y, dim, normalized);
%  CrossCorr = crosscorrfft(x, y, dim, normalized, vector);
%
%  Input:
%     x              array of data
%     y              array of data
%     dim            integer
%                      Dimension over which to perform FFT, first
%                      non-singleton dimension by default
%     normalized     1: normalize by the variance (default)
%                    0: no normalization
%     vector         integer
%                      Dimension over which to sum as components of a
%                      vector (scalar behavior by default)
%
%  Output:
%     CrossCorr       autocorrelation of y

function CrossCorr = crosscorrfft(varargin)

switch nargin
  case 1
    % normalized output by default
    x = varargin{1};
    y = varargin{2};
    dim = [];
    normalized = 1;
    vector = 0;
  case 2
    x = varargin{1};
    y = varargin{2};
    dim = varargin{3};
    normalized = 1;
    vector = 0;
  case 3
    x = varargin{1};
    y = varargin{2};
    dim = varargin{3};
    normalized = varargin{4};
    vector = 0;
  case 4
    x = varargin{1};
    y = varargin{2};
    dim = varargin{3};
    normalized = varargin{4};
    vector = varargin{5};
  otherwise
    error('Wrong number of input arguments.')
end

Dimsx = ndims(x);
Dimsy = ndims(y);

if ~isequal(Dimsx,Dimsy)
  error('Input arrays must be the same size.')
end

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
x = bsxfun(@minus, x, mean(x, dim));
y = bsxfun(@minus, y, mean(y, dim));


N = size(y, dim);
Fx = fft(x, 2*N, dim);
Fy = fft(y, 2*N, dim);
r = ifft(Fx.*conj(Fy), [], dim);

% select only the first half of the FFT dimension
idx = cell(1, Dimsy);
idx(:) = {':'};
idx(dim) = {1:N};
CrossCorr = real(r(idx{:}));

if vector
  % sum along dimension given by value of vector argument
  CrossCorr = sum(CrossCorr, vector);
end

n = N*ones(1, N) - [1:N] + 1;
sizen = ones(1,Dimsy);
sizen(dim) = N;
n = reshape(n,sizen);
CrossCorr = CrossCorr./n;

if normalized
  % divide by first value of the FFT dimension
  idx = cell(1, Dimsy);
  idx(:) = {':'};
  idx{dim} = 1;
  CrossCorr = CrossCorr./CrossCorr(idx{:});
end

end