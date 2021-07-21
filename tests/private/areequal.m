function ok = areequal(A,B,threshold,mode)

if ~isnumeric(A) || ~isnumeric(B)
  error('Both inputs to areequal() must be numeric.');
end

if numel(A)~=numel(B) || ndims(A)~=ndims(B) || any(size(A)~=size(B))
    error('The two inputs to areequal() must be arrays of the same size.');
end

if nargin<3
  threshold = 0;
end

if nargin<4
  if threshold==0
    mode = 'abs';
  else
    error('Mode (4th input) required, must be either ''rel'' (relative) or ''abs'' (absolute).');
  end
end

if ischar(threshold)
  error('Threshold (3rd input) must be a scalar.');
end

if ~ischar(mode)
  error('Mode (4th input) must be ''abs'' or ''rel''.');
end

delta = abs(A(:)-B(:));
switch mode
  case 'abs'
    ok = all(delta<=threshold);
  case 'rel'
    maxB = max(abs(B(:)));
    ok = all(delta<=threshold*maxB);
  otherwise
    error('Unreognized comparison mode %s',mode);
end

end
