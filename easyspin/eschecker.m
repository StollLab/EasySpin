% ESCHECKER   EasySpin expiry checker
%
%    eschecker
%
% Checks whether the installed version of EasySpin has expired.

% To utilize this expiry checker, add the following lines between %{ and %} to
% a file that should check the expiry date:
%{
% Check expiry date
error(eschecker);
%}

function varargout = eschecker(varargin)

if nargout>1
  error('At most one output argument is possible.');
end

% Determine whether this function has been called by an EasySpin function
EasySpinFilenames = {'cardamom','chili','curry','esfit','garlic','pepper','saffron','salt','spidyan'};
db = dbstack;
internalCall = numel(db)>=2 && any(strcmp(db(2).name,EasySpinFilenames));


% Version and expiry date settings
%-------------------------------------------------------------------------------
persistent esVersion expiryDate
if isempty(esVersion)
  info = easyspin_info;
  esVersion = info.Version;
  expiryDate = info.ExpiryDate;
end


% Don't check at every iteration
%-------------------------------------------------------------------------------
persistent versionOK nCalls
if isempty(versionOK), versionOK = false; end
if isempty(nCalls), nCalls = 0; end
if internalCall
  if versionOK && nCalls<100
    nCalls = nCalls + 1;
    if nargout==1
      varargout = {''};
    else
      varargout = {};
    end
    return
  else
    nCalls = 1;
  end
end


% Check for expiry
%-------------------------------------------------------------------------------
isUnpackaged = esVersion(1)=='$';
if isUnpackaged
  daysToExpiry = Inf;
  versionOK = true;
else
  nowDate = datetime('today');
  expiryDate = datetime(expiryDate);
  daysToExpiry = days(expiryDate-nowDate);
  versionOK = daysToExpiry>0;
end


% Set error message if version is expired
%-------------------------------------------------------------------------------
err = '';
if ~isUnpackaged
  if ~versionOK
    err = sprintf('Your EasySpin version %s has expired.\nPlease visit easyspin.org and update to the latest version.',esVersion);
  end
end


% Show warning if close to expiry date
%-------------------------------------------------------------------------------
if ~isUnpackaged
  if versionOK && daysToExpiry<30
    disp(   '***************************************************************');
    fprintf('This EasySpin version will expire in %d days.\n',daysToExpiry);
    disp(   'Please visit easyspin.org and download the latest version.');
    disp(   '***************************************************************');
  end
end


% Assign or display
%-------------------------------------------------------------------------------
if internalCall
  varargout = {err};
else
  if nargout==1
    varargout = {err};
  else
    if isUnpackaged
      fprintf('Your EasySpin installation is valid.\nIt is unpackaged.\n');
    elseif versionOK
      fprintf('Your EasySpin installation is valid.\nIt expires in %d days.\n',daysToExpiry);
    else
      disp(err);
    end
  end
end
