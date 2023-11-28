function ok = test()

% Test whether esfit_legacy() returns correct confidence intervals (CIs).
% This test compares esfit_legacy()'s CI values to CI values obtained from
% MATLAB's cftool and DeerLab's fitparamodel (which are consistent).

data = textread('esfit_std_data.txt'); %#ok<DTXTRD>
t = data(:,1);
y = data(:,2);

pmodel = [2 0.3 0.2];
lb = [0.1 0.1 -2];
p0 = [1.5 0.3 0.5];
ub = [10 0.5 2];

modelfcn = @(p) p(1)*exp(-p(2)*t)+p(3);

FitOpt = struct;
FitOpt.Method = 'simplex fcn';
FitOpt.Verbosity = 0;
FitOpt.AutoScale = 'none';

result = esfit_legacy(y,modelfcn,p0,lb,ub,FitOpt);

% Reference parameter standard deviations, determined using MATLAB's
% Curve Fitting Tool (cftool) and DeerLab's fitparamodel.
pstd_ref = [0.0149, 0.0043, 0.0053];

ok(1) = areequal(result.pfit,pmodel.',3e-2,'abs');
ok(2) = areequal(result.pstd,pstd_ref.',2e-2,'rel');

end
