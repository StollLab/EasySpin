% apowin  Apodization windows 
%
%    w = apowin(Type,n)
%    w = apowin(Type,n,alpha)
%
%    Returns an apodization window. n is the number
%    of points. The string Type specifies the type
%    and can be
%
%      'bla'    Blackman
%      'bar'    Bartlett
%      'con'    Connes
%      'cos'    Cosine
%      'ham'    Hamming
%      'han'    Hann (also called Hanning)
%      'wel'    Welch
%
%   The following three windows need the parameter
%   alpha. Reasonable ranges for alpha are given.
%
%      'exp'    Exponential    2 to 6
%      'gau'    Gaussian       0.6 to 1.2
%      'kai'    Kaiser         3 to 9
%
%   A '+' ('-') appended to Type indicates that only the
%   right (left) half of the window should be constructed.
%
%      'ham'    symmetric (-1 <= x <= 1, n points)
%      'ham+'   right side only (0 <= x <= 1, n points)
%      'ham-'   left side only (-1 <= x <= 0, n points)

function varargout = apowin(Type,n,alpha)

if (nargin==0), help(mfilename); return; end

if ~ischar(Type) | numel(Type)<3 | numel(Type)>4
  error('The first argument must be a 3- or 4-character string!');
end

if numel(Type)==4,
  switch Type(4)
    case '+', xmin = 0; xmax = 1;
    case '-', xmin = -1; xmax = 0;
    otherwise
      error('Wrong 4th character in Type. Should be + or -.');
  end
else
  xmin = -1; xmax = 1;
end

if ~isreal(n) | ~(n>0) | mod(n,1)
  error('n must be a positive integer!');
end

x = linspace(xmin,xmax,n);

% Window type switchyard
%-------------------------------------------------------
Type = Type(1:3);
n_arg = 2;
if strcmp(Type,'ham')
  n_arg = 2;
  w = 0.54 + 0.46*cos(pi*x);
  
elseif strcmp(Type,'kai')
  n_arg = 3;
  w = besseli(0,alpha*sqrt(1-x.^2))/besseli(0,alpha);
  
elseif strcmp(Type,'gau')
  n_arg = 3;
  w = exp(-2*x.^2/alpha^2);
  
elseif strcmp(Type,'exp')
  n_arg = 3;
  w = exp(-alpha*abs(x));
  
elseif strcmp(Type,'han')
  n_arg = 2;
  w = 0.5 + 0.5*cos(pi*x);
  
elseif strcmp(Type,'bla')
  n_arg = 2;
  w = 0.42 + 0.5*cos(pi*x) + 0.08*cos(2*pi*x);
  
elseif strcmp(Type,'bar')
  n_arg = 2;
  w = 1 - abs(x);
  
elseif strcmp(Type,'con')
  n_arg = 2;
  w = (1-x.^2).^2;
  
elseif strcmp(Type,'cos')
  n_arg = 2;
  w = cos(pi*x/2);
  
elseif strcmp(Type,'wel')
  n_arg = 2;
  w = 1 - x.^2;
  
else
  error('Unknown apodization window specified!');
  
end
if (nargin~=n_arg)
  error('Wrong number of input arguments');
end

if (xmax-xmin==2),
  % Symmetrize (remove possible numerical asymmetries)
  w = (w+w(end:-1:1))/2;
end

w = w/max(w); % normalize maximum to 1

w = w(:); % convert to column vector

switch nargout
  case 0
    plot(w);
  case 1
    varargout = {w};
end

return
