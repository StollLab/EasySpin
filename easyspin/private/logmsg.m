% logmsg    EasySpin logging 
%
%   logmsg(MsgLevel,varargin)
%
%   Prints log messages to the screen, if the global EasySpinLogLevel
%   allows. varargin is passed on to fprintf.
%
%   Log levels:
%     0  no messages
%     1  normal messages
%     2  detailed messages
%     3  loop counters
%     4  messages from inside loops (debug level)
%
%   A message is only displayed if MsgLevel is smaller or equal to 
%   EasySpinLogLevel.

function logmsg(varargin)

% Connect to global variable
global EasySpinLogLevel

if isempty(EasySpinLogLevel), return; end

MsgLevel = varargin{1};

if numel(MsgLevel)~=1 || ~isnumeric(MsgLevel) || mod(MsgLevel,1)~=0 || MsgLevel<0
  error('The first input to logmsg must be a non-negative integer.');
end

% Don't display if message level is above loglevel
if MsgLevel>EasySpinLogLevel
  return
end

% Display message
if nargin>1
  fprintf(varargin{2:end});
  fprintf('\n');
else
  error('At least two input arguments expected!');
end
