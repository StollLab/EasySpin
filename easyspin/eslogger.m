% logmsg    EasySpin logging function
%
%   eslogger(logLevel)             % set log level
%   logLevel = logmsg            % get log level
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
%     4  (trace) very detailed debug messages

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
    return
  else
    error('The log level must be 0, 1, 2, 3 or 4.');
  end
end

msgLevel = varargin{1};
args = varargin(2:end);

if ~isValidLogLevel(msgLevel)
  error('The message log level must be 0, 1, 2, 3 or 4.');
end

% Don't display if message level is above loglevel
if msgLevel>logLevel
  return
end

% Display message
if numel(args)>=1
  fprintf(args{:});
  fprintf('\n');
end

end

function tf = isValidLogLevel(lev)
tf = isnumeric(lev) && isscalar(lev) && ismember(lev,0:4);
end
