% ESCHECKER   EasySpin expiry checker
%
%    eslicense
%
%   Checks the expiry of EasySpin

function varargout = eschecker(varargin)

if (nargout>0), varargout = cell(1,nargout); end

% Add the following lines to a file that should check the license:
% --------License ------------------------------------------------
%
%LicErr = 'Could not determine license.';
%Link = 'epr@eth'; eschecker; error(LicErr); clear Link LicErr;
%
% --------License ------------------------------------------------
Link1 = 'epr@eth';
Link2 = 'LicErr';

try
  Link = evalin('caller','Link');
catch
  Link = '';
end
InternalCall = strcmp(Link,Link1);

persistent LicenseOK nCalls;
if isempty(LicenseOK), LicenseOK = 0; end
if isempty(nCalls), nCalls = 1; end

if InternalCall
  if LicenseOK && (nCalls<100)
    nCalls = nCalls + 1;
    assignin('caller',Link2,'');
    return
  else
    nCalls = 1;
  end
end

% Time-limited license
%------------------------------------------------------------------
% License not valid between ExpiryDate
% file in tempdir with dates beyond HorizonDate are ignored
ExpiryDate = 888888;
HorizonDate = 999999;

% Determine modification date of temp directory
%--------------------------------------------------------------------
try
  List = dir(tempdir); % List(1).name should be '.'
  if numel(List(1).name)==1
    TmpDirDate = datenum(List(1).date);
  else
    TmpDirDate = [];
  end
catch
  lasterr('');
  TmpDirDate = [];
end
if (TmpDirDate==0), TmpDirDate = []; end
if (TmpDirDate>=HorizonDate), TmpDirDate = []; end


% Determine system time
%--------------------------------------------------------------------
NowDate = datenum(builtin('clock'));


% Expired or not?
%--------------------------------------------------------------------
Expired = (NowDate>ExpiryDate);
if ~isempty(TmpDirDate)
  Expired = Expired | any(TmpDirDate>ExpiryDate);
end

if (Expired)
  err = 'Your EasySpin version has expired. Please visit easyspin.org and update to the latest version.';
else
  DaysToExpiry = fix(ExpiryDate - NowDate);
  if (DaysToExpiry<30)
    disp(   '***************************************************************');
    fprintf('This EasySpin version will expire in %d days.\n',DaysToExpiry);
    disp(   'Please visit easyspin.org and download the latest version.');
    disp(   '***************************************************************');
  end
  err = '';
end
LicenseOK = isempty(err);

% assign or display
%--------------------------------------------------------------------
if InternalCall

  assignin('caller',Link2,err);

else

  if LicenseOK
    disp('Your EasySpin installation is ok.');
  else
    disp(err);
  end
end
