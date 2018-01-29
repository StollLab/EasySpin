function varargout = runprivate2(fname,varargin)

result = cell(1,nargout);

[result{:}] = feval(fname,varargin{:});

varargout = result;

return
