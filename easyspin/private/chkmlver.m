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
  ErrorString = sprintf('\n  ================================================================\n  Easyspin: MATLAB version not supported\n  ----------------------------------------------------------------\n  This version of EasySpin requires MATLAB %s and later.\n  Your MATLAB version %s is not supported.\n  ================================================================\n',MinimalVersionString,version);
else
  ErrorString = '';
end
