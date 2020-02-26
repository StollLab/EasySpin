function ok = test()
%{
The following information is contained in the versions.txt file that is 
found in the data subfolder of the testing directory and called by
easyspinversioncheck if called with Opt.Test = true:

stable:5.2.25
default:7.1.2-alpha.3
development:6.0.0-dev.8

A total of 14 tests is run.

%}

Opt.Silent = true;
Opt.Test = true;

Errors = zeros(11,1);

% -------------------------------------------------------------
% Test 1:
VInfo.Version = '5.2.2';
VInfo.ReleaseChannel = 'stable';

[TF, Version] = easyspinversioncheck(VInfo,Opt);

if TF == 0 || ~strcmp(Version,'5.2.25')
  Errors(1) = true;
end

% -------------------------------------------------------------
% Test 2:
VInfo.Version = '5.2.2';
VInfo.ReleaseChannel = 'development';

[TF, Version] = easyspinversioncheck(VInfo,Opt);

if TF == 0 || ~strcmp(Version,'6.0.0-dev.8')
  Errors(2) = true;
end

% -------------------------------------------------------------
% Test 3:
VInfo.Version = '5.2.2';
VInfo.ReleaseChannel = 'not a branch';

[TF, Version] = easyspinversioncheck(VInfo,Opt);

if TF == 1 || ~isempty(Version)
  Errors(3) = true;
end

% -------------------------------------------------------------
% Test 4:
VInfo.Version = '5.3.2';
VInfo.ReleaseChannel = 'stable';

[TF, Version] = easyspinversioncheck(VInfo,Opt);

if TF == 1 || ~strcmp(Version,'5.2.25')
  Errors(4) = true;
end

% -------------------------------------------------------------
% Test 5:
VInfo.Version = '6.3.2';
VInfo.ReleaseChannel = 'stable';

[TF, Version] = easyspinversioncheck(VInfo,Opt);

if TF == 1 || ~strcmp(Version,'5.2.25')
  Errors(5) = true;
end

% -------------------------------------------------------------
% Test 5:
VInfo.Version = '6.3.2';
VInfo.ReleaseChannel = 'development';

[TF, Version] = easyspinversioncheck(VInfo,Opt);

if TF == 1 || ~strcmp(Version,'6.0.0-dev.8')
  Errors(6) = true;
end

% -------------------------------------------------------------
% Test 7:
VInfo.Version = '6.3.2-alpha.1';
VInfo.ReleaseChannel = 'development';

[TF, Version] = easyspinversioncheck(VInfo,Opt);

if TF == 1 || ~strcmp(Version,'6.0.0-dev.8')
  Errors(7) = true;
end

% -------------------------------------------------------------
% Test 8:
VInfo.Version = '6.0.0';
VInfo.ReleaseChannel = 'development';

[TF, Version] = easyspinversioncheck(VInfo,Opt);

if TF == 1 || ~strcmp(Version,'6.0.0-dev.8')
  Errors(8) = true;
end

% -------------------------------------------------------------
% Test 9:
VInfo.Version = '6.0.0-dev.1';
VInfo.ReleaseChannel = 'development';

[TF, Version] = easyspinversioncheck(VInfo,Opt);

if TF == 0 || ~strcmp(Version,'6.0.0-dev.8')
  Errors(9) = true;
end

% -------------------------------------------------------------
% Test 10:
VInfo.Version = '6.0.0-alpha.1';
VInfo.ReleaseChannel = 'development';

[TF, Version] = easyspinversioncheck(VInfo,Opt);

if TF == 1 || ~strcmp(Version,'6.0.0-dev.8')
  Errors(10) = true;
end


% -------------------------------------------------------------
% Test 11:
VInfo.Version = '$ReleaseID$';
VInfo.ReleaseChannel = '$ReleaseChannel$';

[TF, Version] = easyspinversioncheck(VInfo,Opt);

if TF == 1 || ~isempty(Version)
  Errors(11) = true;
end

% -------------------------------------------------------------
% Test 12:
VInfo.Version = '7.1.2-dev.10';
VInfo.ReleaseChannel = 'default';

[TF, Version] = easyspinversioncheck(VInfo,Opt);

if TF == 0 || ~strcmp(Version,'7.1.2-alpha.3')
  Errors(12) = true;
end

% -------------------------------------------------------------
% Test 13:
VInfo.Version = '7.1.2-beta';
VInfo.ReleaseChannel = 'default';

[TF, Version] = easyspinversioncheck(VInfo,Opt);

if TF == 1 || ~strcmp(Version,'7.1.2-alpha.3')
  Errors(13) = true;
end

% -------------------------------------------------------------
% Test 14:
VInfo.Version = '7.1.2-alpha.12';
VInfo.ReleaseChannel = 'default';

[TF, Version] = easyspinversioncheck(VInfo,Opt);

if TF == 1 || ~strcmp(Version,'7.1.2-alpha.3')
  Errors(14) = true;
end

% -------------------------------------------------------------
% Validation

ok = ~any(Errors);
