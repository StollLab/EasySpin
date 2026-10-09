function ok = test()

% Strain widths from resfields_perturb against widths computed from finite
% differences of the perturbation-theory resonance fields

[S,c] = perturbstrainsystems();

Exp.mwFreq = 9.5;
Exp.Range = [100 1000];
Exp.SampleFrame = [0.2 0.9 0.4; 1.3 0.4 2.1];

ok = true(2,c);
for order = 1:2
  Opt.PerturbOrder = order;
  for k = 1:c
    Sys = S{k};
    [~,~,W] = resfields_perturb(Sys,Exp,Opt);
    Sys_ = runprivate('validatespinsys',Sys);
    Cov = Sys_.StrainData.Q*Sys_.StrainData.Q.';
    pars = Sys.StrainPars;
    Sys = rmfield(Sys,intersect(fieldnames(Sys),{'StrainPars','StrainFWHM','StrainCorr','StrainModes'}));

    % HStrain contribution, added in quadrature
    wH2 = 0;
    if isfield(Sys,'HStrain')
      [~,~,WH] = resfields_perturb(Sys,Exp,Opt);
      wH2 = WH.^2;
    end

    J = zeros([size(W) numel(pars)]);
    for p = 1:numel(pars)
      [f,idx] = parseref(pars{p});
      h = 1e-6*max(1,abs(Sys.(f)(idx)));
      Sp = Sys; Sp.(f)(idx) = Sp.(f)(idx) + h;
      Sm = Sys; Sm.(f)(idx) = Sm.(f)(idx) - h;
      J(:,:,p) = (resfields_perturb(Sp,Exp,Opt) - resfields_perturb(Sm,Exp,Opt))/(2*h);
    end
    J = reshape(J,[],numel(pars));
    Wfd = reshape(sum((J*Cov).*J,2),size(W));
    Wfd = sqrt(Wfd + wH2);
    if isempty(W), W = zeros(size(Wfd)); end
    ok(order,k) = areequal(W,Wfd,1e-5,'rel');
  end
end
ok = ok(:).';

% Widths without the frequency-to-field conversion, in MHz
Sys = S{2};
[~,~,W1] = resfields_perturb(Sys,Exp);
[~,~,WH1] = resfields_perturb(struct('g',Sys.g,'gFrame',Sys.gFrame,'HStrain',[1 1 1]),Exp);
Opt = struct('Freq2Field',0);
[~,~,W0] = resfields_perturb(Sys,Exp,Opt);
[~,~,WH0] = resfields_perturb(struct('g',Sys.g,'gFrame',Sys.gFrame,'HStrain',[1 1 1]),Exp,Opt);
ok(end+1) = areequal(W0,W1./(WH1(1,:)./WH0(1,:)),1e-10,'rel');

% Negative fields: mirrored resonances have the same widths
Exp.Range = [-2000 2000];
[B,~,W] = resfields_perturb(S{2},Exp);
n = size(B,1)/2;
ok(end+1) = areequal(B(n+1:end,:),-B(1:n,:),1e-10,'rel') && ...
  areequal(W(n+1:end,:),W(1:n,:),1e-10,'rel');

end

%-------------------------------------------------------------------------------
function [f,idx] = parseref(str)
tok = regexp(str,'^(\w+)','tokens','once');
f = tok{1};
idx = sscanf(str(numel(f)+1:end),'(%d)');
if isempty(idx), idx = 1; end
end
