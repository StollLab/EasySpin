%  autocorrfft  Calculate autocorrelation function using the FFT.
%
%  AutoCorr = autocorrfft(x);
%  AutoCorr = autocorrfft(x, dim);
%  AutoCorr = autocorrfft(x, dim, normalized);
%  AutoCorr = autocorrfft(x, dim, normalized, vector);
%  AutoCorr = autocorrfft(x, dim, normalized, vector, centered);
%
%  Input:
%     x              array of data
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
%     AutoCorr       autocorrelation of x

function AutoCorr = autocorrfft(varargin)

switch nargin
  case 1
    % normalized output by default
    x = varargin{1};
    dim = [];
    normalized = true;
    vector = 0;
    centered = true;
  case 2
    x = varargin{1};
    dim = varargin{2};
    normalized = true;
    vector = 0;
    centered = true;
  case 3
    x = varargin{1};
    dim = varargin{2};
    normalized = varargin{3};
    vector = 0;
    centered = true;
  case 4
    x = varargin{1};
    dim = varargin{2};
    normalized = varargin{3};
    vector = varargin{4};
    centered = true;
  case 5
    x = varargin{1};
    dim = varargin{2};
    normalized = varargin{3};
    vector = varargin{4};
    centered = varargin{5};
  otherwise
    error('Wrong number of input arguments.')
end

AutoCorr = crosscorrfft(x, x, dim, normalized, vector, centered);

end