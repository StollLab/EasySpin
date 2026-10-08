function ok = test()

% Covariance input for strains: StrainFWHM/StrainCorr and StrainModes

Sys.S = 1;
Sys.g = [2 2.1 2.2];
Sys.D = [300 50];
Sys.StrainPars = {'g(1)','g(2)','D(1)','D(2)'};
fwhm = [0.01 0.02 30 10];
C = [1 0.5 -0.2 0.1; 0.5 1 0 0.3; -0.2 0 1 0.4; 0.1 0.3 0.4 1];
Cov = diag(fwhm)*C*diag(fwhm);

% Full and upper-triangular correlation matrix (row order 12,13,14,23,24,34)
Sys.StrainFWHM = fwhm;
Sys.StrainCorr = C;
S1 = runprivate('validatespinsys',Sys);
Sys.StrainCorr = [0.5 -0.2 0.1 0 0.3 0.4];
S2 = runprivate('validatespinsys',Sys);
ok(1) = areequal(S1.StrainData.Q*S1.StrainData.Q.',Cov,1e-10,'rel');
ok(2) = areequal(S2.StrainData.Q*S2.StrainData.Q.',Cov,1e-10,'rel');

% Mode vectors, including a sign flip of one vector
Sys = rmfield(Sys,{'StrainFWHM','StrainCorr'});
Q = S1.StrainData.Q;
Sys.StrainModes = Q.';
S3 = runprivate('validatespinsys',Sys);
Sys.StrainModes(1,:) = -Sys.StrainModes(1,:);
S4 = runprivate('validatespinsys',Sys);
ok(3) = areequal(S3.StrainData.Q*S3.StrainData.Q.',Cov,1e-10,'rel');
ok(4) = areequal(S4.StrainData.Q*S4.StrainData.Q.',Cov,1e-10,'rel');

% Default: uncorrelated
Sys = rmfield(Sys,'StrainModes');
Sys.StrainFWHM = fwhm;
S5 = runprivate('validatespinsys',Sys);
ok(5) = areequal(S5.StrainData.Q*S5.StrainData.Q.',diag(fwhm.^2),1e-10,'rel');

% Errors
badCorr = {[1 0.5 0; 0.4 1 0; 0 0 1 ], ... % wrong size
           [0.5 2 0 0 0 0], ...             % |c|>1
           [0.9 0.9 -0.9 0 0 0], ...        % not positive semidefinite
           2*eye(4), ...                    % diagonal not 1
           [], ...                          % empty
           [1 0.5 0 0; 0.4 1 0 0; 0 0 1 0; 0 0 0 1]}; % not symmetric
for k = 1:numel(badCorr)
  Sys.StrainCorr = badCorr{k};
  [~,err] = runprivate('validatespinsys',Sys);
  ok(end+1) = ~isempty(err); %#ok<*AGROW>
end
Sys = rmfield(Sys,'StrainCorr');

Sys.StrainFWHM = [1 2 3];
[~,err] = runprivate('validatespinsys',Sys); ok(end+1) = ~isempty(err);
Sys.StrainFWHM = [1 2 3 -4];
[~,err] = runprivate('validatespinsys',Sys); ok(end+1) = ~isempty(err);
Sys.StrainFWHM = fwhm;
Sys.StrainModes = Q.';
[~,err] = runprivate('validatespinsys',Sys); ok(end+1) = contains(err,'either');
Sys = rmfield(Sys,'StrainFWHM');
Sys.StrainCorr = C;
[~,err] = runprivate('validatespinsys',Sys); ok(end+1) = ~isempty(err);

% n = 1: StrainCorr must be 1, empty is an error
Sys = struct('g',[2 2 2.2],'StrainPars',{{'g(3)'}},'StrainFWHM',0.01,'StrainCorr',1);
[~,err] = runprivate('validatespinsys',Sys); ok(end+1) = isempty(err);
Sys.StrainCorr = [];
[~,err] = runprivate('validatespinsys',Sys); ok(end+1) = contains(err,'empty');
