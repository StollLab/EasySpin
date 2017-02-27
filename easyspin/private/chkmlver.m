% chkmlver  Check MATLAB version for EasySpin
%
%    ErrStr = chkmlver
%
%    Returns a non-empty error string in ErrStr if the Matlab
%    version is not suitable for running EasySpin.

function ErrorString = chkmlver

MinimalVersion = '7.5';
MinimalVersionString = '7.5 (R2007b)';

if verLessThan('matlab',MinimalVersion)
  ErrorString = sprintf('\n  Easyspin problem\n  =======================================================\n  EasySpin supports Matlab %s and later.\n  Your Matlab version %s is too old.\n  Please install a more recent version.',MinimalVersionString,version);
else
  ErrorString = '';
end
