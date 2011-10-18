function varargout = runprivate(fname,varargin)

q = feval(fname,varargin{:});

if iscell(q)
  varargout = q;
else
  varargout = {q};
end

return
