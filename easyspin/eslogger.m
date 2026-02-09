% eslogger    EasySpin logging function
%
%   eslogger(logLevel)             % set log level
%   logLevel = eslogger            % get log level
%
%   eslogger(msgLevel,varargin)    % log message
%
%   Prints log messages to the command window if the log level of the message
%   is lower or equal to the set log level. varargin is passed on to fprintf.
%
%   Log levels:
%     0  (off) no messages
%     1  (info) normal messages
%     2  (detail) detailed messages
%     3  (debug) debug messages

function varargout = eslogger(varargin)

% Store log level in persistent variable
persistent logLevel
if isempty(logLevel)
  logLevel = 0;  % default
end

% Get log level
if nargin==0
  varargout = {logLevel};
  return
end

% Set log level
if nargin==1
  newLogLevel = varargin{1};
  if isValidLogLevel(newLogLevel)
    logLevel = newLogLevel;
  else
    error('The log level must be 0, 1, 2, or 3.');
  end
  if logLevel>=1
    printMessage('EasySpin log level set to %d.',logLevel);
  end
  return
end

msgLevel = varargin{1};
args = varargin(2:end);

if ~isValidLogLevel(msgLevel)
  error('The message log level must be 0, 1, 2, or 3.');
end

% Don't log if message level is above log level
if msgLevel>logLevel
  return
end

printMessage(args{:});

end


function tf = isValidLogLevel(lev)
tf = isnumeric(lev) && isscalar(lev) && ismember(lev,[0 1 2 3]);
end


function printMessage(varargin)

% Print function name
printFunctionName = true;
if printFunctionName
  db = dbstack;
  if numel(db)>=3
    funcname = db(3).name;
    lineno = db(3).line;
  else
    funcname = '...';
    lineno = '';
  end
  fprintf('[%s:%d] ',funcname,lineno);
end

% Print message
fprintf(varargin{:});
fprintf('\n');

end
