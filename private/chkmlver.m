% chkmlver  Check MATLAB version for EasySpin
%
%    ErrStr = chkmlver
%
%    Returns a non-empty error string in ErrStr if the Matlab
%    version is not suitable for running EasySpin.

function ErrorString = chkmlver

ErrorString = '';

% From Matlab 5.3.1 to 6, eig() changed its behaviour for
% hermitian matrices. From 6 on, eigenvals are real and sorted,
% before they are complex sometimes and unsorted.
% In Matlab 7.0.0, sort() had to be introduced again, because
% Matlab 7.0.0 did't return eigenvalues in sorted order for
% some quasidiagonal matrices. This Matlab bug was fixed a
% few months later in 7.0.1.

%MinimalVersionString = '6.5 (R13)'; % up to EasySpin 3.1.7
MinimalVersionString = '7.0 (R14)';
MinimalVersion = 7.0;

% Extracts major+minor version number, e.g. 6.0, 6.1, 6.5, 7.0,
% 7.1, 7.2, 7.3, etc. Third and fourth number from version
% ID are not used. E.g. 7.0.0.19901 (R14)  ->  7.0
% 7.12.0.635 (R2011a) -> 7.12
MatlabVersion = sscanf(version,'%f');

if isempty(MatlabVersion)
  ErrorString = 'Could not determine Matlab version.';
  return
end

MatlabVersion = MatlabVersion(1);

if (MatlabVersion<MinimalVersion)
  ErrorString = sprintf('\n  Easyspin problem\n  =======================================================\n  EasySpin supports Matlab %s and later.\n  Your Matlab version %s is too old.\n  Please install a more recent version.',MinimalVersionString,version);
  return
end

% Matlab 7.11 (R2010b) complains about p-file generated with Matlab
% prior to 7.5 (R2007b). As long as we support <R2007b, this warning
% has to switched off.
if MatlabVersion>=7.11
  warning('off','MATLAB:oldPfileVersion');
end
