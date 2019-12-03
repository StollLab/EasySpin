% ESCHECKER   EasySpin expiry checker
%
%    eschecker
%
% Checks whether the installed version of EasySpin has expired.

% Add the following lines between %{ and %} to a file that should check the
% expiry date:
%{
% Check expiry date
error(eschecker);
%}

function varargout = eschecker(varargin)

if nargout>1
  error('At most one output argument is possible.');
end

% Determine whether this function has been called by an EasySpin function
ESfilenames = {'cardamom','chili','curry','esfit','garlic','pepper','saffron','salt','spidyan'};
db = dbstack;
internalCall = numel(db)>=2 && any(strcmp(db(2).name,ESfilenames));

persistent VersionOK nCalls
if isempty(VersionOK), VersionOK = false; end
if isempty(nCalls), nCalls = 1; end

if internalCall
  if VersionOK && nCalls<100
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


% Expiry time settings
%-------------------------------------------------------------------------------
esVersion = '$ReleaseID$';
% Installed version is not valid after ExpiryDate.
ExpiryDate = 888888; % serial date number of expiry date
% Files in tempdir with dates beyond HorizonDate are ignored during the check
HorizonDate = 999999; % serial date number


% Test 1: Expired if system date is after expiry date
%-------------------------------------------------------------------------------
NowDate = fix(datenum(builtin('clock')));
DaysToExpiry = fix(ExpiryDate - NowDate);
isExpired1 = DaysToExpiry<0;


% Test 2: Expired if temp directory modification date is after expiry date
%-------------------------------------------------------------------------------
try
  List = dir(tempdir); % List(1).name should be '.'
  if numel(List(1).name)==1
    TmpDirDate = fix(List(1).datenum);
  else
    TmpDirDate = 0;
  end
catch
  lasterr('');
  TmpDirDate = 0;
end
if TmpDirDate>=HorizonDate, TmpDirDate = 0; end
isExpired2 = TmpDirDate>ExpiryDate;

VersionOK = ~isExpired1 && ~isExpired2;


% Show warning if close to expiry date
%-------------------------------------------------------------------------------
if VersionOK && DaysToExpiry<30
  disp(   '***************************************************************');
  fprintf('This EasySpin version will expire in %d days.\n',DaysToExpiry);
  disp(   'Please visit easyspin.org and download the latest version.');
  disp(   '***************************************************************');
end


% Set error message if version is expired
%-------------------------------------------------------------------------------
if VersionOK
  err = '';
else
  if isExpired1
    err = sprintf('Your EasySpin version %s has expired (1;%d,%d).\nPlease visit easyspin.org and update to the latest version.',esVersion,NowDate,ExpiryDate);
  elseif isExpired2
    err = sprintf('Your EasySpin version %s has expired (2;%d,%d).\nPlease visit easyspin.org and update to the latest version.',esVersion,TmpDirDate,ExpiryDate);
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
    if VersionOK
      fprintf('Your EasySpin installation is ok. It will expire in %d days.\n',DaysToExpiry);
    else
      disp(err);
    end
  end
end
