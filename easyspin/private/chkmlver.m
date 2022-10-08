% chkmlver  Check MATLAB version for EasySpin
%
%    ErrorString = chkmlver
%
%    Returns a non-empty error string in ErrorString if the MATLAB
%    version is not suitable for running EasySpin.
%
%    Usage: error(chkmlver)

function ErrorString = chkmlver

MinimalVersion = '9.2';
MinimalVersionString = '9.2 (R2017a)';

if verLessThan('matlab',MinimalVersion)
  ErrorString = sprintf('\n  Easyspin MATLAB version support\n  =======================================================\n  EasySpin supports MATLAB %s and later.\n  Your MATLAB version %s is not supported.\n.',MinimalVersionString,version);
else
  ErrorString = '';
end
