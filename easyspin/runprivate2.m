function varargout = runprivate2(fname,varargin)

[q1, q2] = feval(fname,varargin{:});

varargout = {q1 q2};

return
