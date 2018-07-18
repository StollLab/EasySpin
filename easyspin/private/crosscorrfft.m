%  crosscorrfft  Calculate cross-correlation function using the FFT.
%
%  CrossCorr = crosscorrfft(x, y);
%  CrossCorr = crosscorrfft(x, y, dim);
%  CrossCorr = crosscorrfft(x, y, dim, normalized);
%  CrossCorr = crosscorrfft(x, y, dim, normalized, vector);
%  CrossCorr = crosscorrfft(x, y, dim, normalized, vector, centered);
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
%     centered       1: subtract mean from data (default)
%                    0: leave data as is
%
%  Output:
%     CrossCorr       autocorrelation of y

function CrossCorr = crosscorrfft(varargin)

switch nargin
  case 2
    % normalized output by default
    x = varargin{1};
    y = varargin{2};
    dim = [];
    normalized = 1;
    vector = 0;
    centered = 1;
  case 3
    x = varargin{1};
    y = varargin{2};
    dim = varargin{3};
    normalized = 1;
    vector = 0;
    centered = 1;
  case 4
    x = varargin{1};
    y = varargin{2};
    dim = varargin{3};
    normalized = varargin{4};
    vector = 0;
    centered = 1;
  case 5
    x = varargin{1};
    y = varargin{2};
    dim = varargin{3};
    normalized = varargin{4};
    vector = varargin{5};
    centered = 1;
  case 6
    x = varargin{1};
    y = varargin{2};
    dim = varargin{3};
    normalized = varargin{4};
    vector = varargin{5};
    centered = varargin{6};
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
if centered
  x = bsxfun(@minus, x, mean(x, dim));
  y = bsxfun(@minus, y, mean(y, dim));
end


N = size(y, dim);
Fx = fft(x, 2*N, dim);
if ~isequal(x,y)
  Fy = fft(y, 2*N, dim);
else
  Fy = Fx;
end
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

n = N*ones(1, N) - (1:N) + 1;
sizen = ones(1,Dimsy);
sizen(dim) = N;
n = reshape(n,sizen);
CrossCorr = bsxfun(@rdivide,CrossCorr,n);

if normalized
  % divide by first value of the FFT dimension
  idx = cell(1, Dimsy);
  idx(:) = {':'};
  idx{dim} = 1;
  CrossCorr = bsxfun(@rdivide,CrossCorr,CrossCorr(idx{:}));
end

end