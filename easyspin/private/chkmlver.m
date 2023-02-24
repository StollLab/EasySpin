% chkmlver  Check MATLAB version for EasySpin
%
%    errorString = chkmlver
%
%    Returns a non-empty error string in errorString if the MATLAB
%    version is not suitable for running EasySpin.
%
%    Usage: error(chkmlver)

function errorString = chkmlver

minimalVersion = '9.2';
minimalVersionString = '9.2 (R2017a)';

if verLessThan('matlab',minimalVersion)
  errorString = sprintf('\n  ================================================================\n  Easyspin: MATLAB version not supported\n  ----------------------------------------------------------------\n  This version of EasySpin requires MATLAB %s and later.\n  Your MATLAB version %s is not supported.\n  ================================================================\n',minimalVersionString,version);
else
  errorString = '';
end
