function varargout = runprivate(fname,varargin)

[q, q2] = feval(fname,varargin{:});

if iscell(q)
  varargout = q;
else
  varargout = {q};
end

return
