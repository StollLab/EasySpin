function varargout = runprivate(fname,varargin)

result = cell(1,nargout);

[result{:}] = feval(fname,varargin{:});

varargout = result;

return
