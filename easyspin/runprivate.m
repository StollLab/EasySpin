% Run a function from the private/ subfolder
%
%   ... = runprivate(fname)
%   ... = runprivate(fname,arg1,arg2,...)
%
% Input:
%   fname       function name (string, character array)
%   arg1, arg2  function arguments
%
% Output:
%   outputs returned from function call

function varargout = runprivate(fname,varargin)

varargout = cell(1,nargout);
[varargout{:}] = feval(fname,varargin{:});

end
