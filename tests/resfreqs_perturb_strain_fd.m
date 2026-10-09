function ok = test()

% Strain widths from resfreqs_perturb against widths computed from finite
% differences of the perturbation-theory resonance frequencies

[S,c] = perturbstrainsystems();

Exp.Field = 350;
Exp.SampleFrame = [0.2 0.9 0.4; 1.3 0.4 2.1];

ok = true(2,c);
for order = 1:2
  Opt.PerturbOrder = order;
  for k = 1:c
    Sys = S{k};
    [~,~,W] = resfreqs_perturb(Sys,Exp,Opt);
    Sys_ = runprivate('validatespinsys',Sys);
    Cov = Sys_.StrainData.Q*Sys_.StrainData.Q.';
    pars = Sys.StrainPars;
    Sys = rmfield(Sys,intersect(fieldnames(Sys),{'StrainPars','StrainFWHM','StrainCorr','StrainModes'}));

    % HStrain contribution, added in quadrature
    wH2 = 0;
    if isfield(Sys,'HStrain')
      [~,~,WH] = resfreqs_perturb(Sys,Exp,Opt);
      wH2 = WH.^2;
    end

    J = zeros([size(W) numel(pars)]);
    for p = 1:numel(pars)
      [f,idx] = parseref(pars{p});
      h = 1e-6*max(1,abs(Sys.(f)(idx)));
      Sp = Sys; Sp.(f)(idx) = Sp.(f)(idx) + h;
      Sm = Sys; Sm.(f)(idx) = Sm.(f)(idx) - h;
      J(:,:,p) = (resfreqs_perturb(Sp,Exp,Opt) - resfreqs_perturb(Sm,Exp,Opt))/(2*h);
    end
    J = reshape(J,[],numel(pars));
    Wfd = reshape(sum((J*Cov).*J,2),size(W));
    Wfd = sqrt(Wfd + wH2);
    if isempty(W), W = zeros(size(Wfd)); end
    ok(order,k) = areequal(W,Wfd,1e-5,'rel');
  end
end
ok = ok(:).';

end

%-------------------------------------------------------------------------------
function [f,idx] = parseref(str)
tok = regexp(str,'^(\w+)','tokens','once');
f = tok{1};
idx = sscanf(str(numel(f)+1:end),'(%d)');
if isempty(idx), idx = 1; end
end
