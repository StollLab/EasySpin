function ok = test()

% Strain widths from resfields (field domain) against widths computed from
% finite differences. For each resonance at field B, the derivatives of the
% resonance field are obtained by implicit differentiation,
%   dB/dp = -(dnu/dp)/(dnu/dB),
% with dnu/dp and dnu/dB from finite differences of resfreqs_matrix at B.

ang = [0.3 0.7 1.1];
c = 0;
c=c+1; S{c} = struct('g',[2 2.1 2.2],'gFrame',ang,'StrainPars',{{'g(1)','g(3)','gFrame(2)'}},'StrainFWHM',[0.01 0.02 0.05],'StrainCorr',[0.6 0 0]);
c=c+1; S{c} = struct('Nucs','63Cu','g',[2.05 2.25],'A',[60 500],'StrainPars',{{'g(2)','A(2)'}},'StrainFWHM',[0.02 30],'StrainCorr',-0.8,'HStrain',[10 10 20]);
c=c+1; S{c} = struct('S',1,'D',[1000 200],'DFrame',ang,'StrainPars',{{'D(1)','D(2)','DFrame(3)'}},'StrainModes',[50 20 0; 0 0 0.05]);
c=c+1; S{c} = struct('S',[1/2 1/2],'g',[2 2.1],'J',100,'dip',20,'StrainPars',{{'g(2)','J'}},'StrainFWHM',[0.01 5]);

Exp.mwFreq = 9.5;
Exp.Range = [100 500];
Exp.SampleFrame = [0.2 0.9 0.4];
Opt.FuzzLevel = 0; % no random noise, for finite differences
Opt.ModellingAccuracy = 1e-12; % accurate eigenvectors at the resonance fields

for k = 1:c
  Sys = S{k};
  [B,~,W,T] = resfields(Sys,Exp,Opt);
  Sys_ = runprivate('validatespinsys',Sys);
  Cov = Sys_.StrainData.Q*Sys_.StrainData.Q.';
  pars = Sys.StrainPars;
  Sys = rmfield(Sys,intersect(fieldnames(Sys),{'StrainPars','StrainFWHM','StrainCorr','StrainModes'}));

  % HStrain contribution, added in quadrature
  wH2 = 0;
  if isfield(Sys,'HStrain')
    [~,~,WH] = resfields(Sys,Exp,Opt);
    wH2 = WH.^2;
  end

  Wfd = zeros(size(W));
  for iRes = 1:numel(B)
    ExpF.SampleFrame = Exp.SampleFrame;
    OptF = Opt;
    OptF.Transitions = T(iRes,:);
    hB = 1e-4;
    ExpF.Field = B(iRes)+hB; nup = resfreqs_matrix(Sys,ExpF,OptF);
    ExpF.Field = B(iRes)-hB; num = resfreqs_matrix(Sys,ExpF,OptF);
    dnudB = (nup-num)/(2*hB);
    ExpF.Field = B(iRes);
    J = zeros(1,numel(pars));
    for p = 1:numel(pars)
      [f,idx] = parseref(pars{p});
      h = 1e-6*max(1,abs(Sys.(f)(idx)));
      Sp = Sys; Sp.(f)(idx) = Sp.(f)(idx) + h;
      Sm = Sys; Sm.(f)(idx) = Sm.(f)(idx) - h;
      J(p) = -(resfreqs_matrix(Sp,ExpF,OptF) - resfreqs_matrix(Sm,ExpF,OptF))/(2*h)/dnudB;
    end
    Wfd(iRes) = J*Cov*J.';
  end
  Wfd = sqrt(Wfd + wH2);
  ok(k) = areequal(W,Wfd,1e-4,'rel');
end

end

%-------------------------------------------------------------------------------
function [f,idx] = parseref(str)
tok = regexp(str,'^(\w+)','tokens','once');
f = tok{1};
idx = sscanf(str(numel(f)+1:end),'(%d)');
if isempty(idx), idx = 1; end
end
