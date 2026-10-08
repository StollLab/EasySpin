function ok = test()

% Strain widths from resfreqs_matrix against widths computed from finite
% differences of the resonance frequencies

ang = [0.3 0.7 1.1];
c = 0;
c=c+1; S{c} = struct('g',[2 2.1 2.2],'gFrame',ang,'StrainPars',{{'g(1)','g(3)','gFrame(2)'}},'StrainFWHM',[0.01 0.02 0.05],'StrainCorr',[0.6 0 0]);
c=c+1; S{c} = struct('g',[2 2.2 2.0 0.01 0.02 0.03],'StrainPars',{{'g(2)','g(6)'}},'StrainFWHM',[0.01 0.005]);
c=c+1; S{c} = struct('Nucs','14N','g',2,'A',[20 50],'AFrame',ang,'Q',[3 0.3],'QFrame',ang,'StrainPars',{{'A(1)','A(2)','AFrame(2)','Q(1)','Q(2)','QFrame(3)'}},'StrainFWHM',[2 5 0.05 0.3 0.1 0.05]);
c=c+1; S{c} = struct('S',1,'D',[3000 500],'DFrame',ang,'StrainPars',{{'D(1)','D(2)','DFrame(1)'}},'StrainFWHM',[100 50 0.05],'StrainCorr',[0.5 0 0],'HStrain',[10 20 30]);
c=c+1; S{c} = struct('S',3/2,'D',[-1000 -500 1500],'StrainPars',{{'D(2)'}},'StrainFWHM',50);
c=c+1; S{c} = struct('S',[1/2 1/2],'g',[2 2.1 2.2; 2.05 2.05 2.0],'J',300,'dip',[20 5],'dvec',[5 0 0],'eeFrame',ang,'StrainPars',{{'J','dip(2)','dvec(1)','eeFrame(2)','g(2,3)'}},'StrainFWHM',[20 2 2 0.05 0.01]);
c=c+1; S{c} = struct('S',[1/2 1/2],'g',[2 2.1],'ee',[200 250 300],'StrainPars',{{'ee(2)'}},'StrainModes',[10; -10]);
c=c+1; S{c} = struct('Nucs','1H,1H','A',[50 80],'nn',[1 2 3],'sigma',[1 1.01 1.02; 1 1 1],'sigmaFrame',[ang; 0 0 0],'StrainPars',{{'A(2)','nn(3)','sigma(1,2)','sigmaFrame(1,2)'}},'StrainFWHM',[3 1 0.01 0.1]);

Exp.Field = 340;
Exp.SampleFrame = [0.2 0.9 0.4];
Opt.Threshold = 0;
Opt.FuzzLevel = 0; % no random noise, for finite differences

for k = 1:c
  Sys = S{k};
  [P,~,W,T] = resfreqs_matrix(Sys,Exp,Opt);
  Opt_ = Opt; Opt_.Transitions = T;
  % Covariance matrix and parameter list
  Sys_ = runprivate('validatespinsys',Sys);
  Cov = Sys_.StrainData.Q*Sys_.StrainData.Q.';
  pars = Sys.StrainPars;
  Sys = rmfield(Sys,intersect(fieldnames(Sys),{'StrainPars','StrainFWHM','StrainCorr','StrainModes'}));
  J = zeros(numel(P),numel(pars));
  for p = 1:numel(pars)
    [f,idx] = parseref(pars{p},Sys);
    h = 1e-6*max(1,abs(Sys.(f)(idx)));
    Sp = Sys; Sp.(f)(idx) = Sp.(f)(idx) + h;
    Sm = Sys; Sm.(f)(idx) = Sm.(f)(idx) - h;
    J(:,p) = (resfreqs_matrix(Sp,Exp,Opt_) - resfreqs_matrix(Sm,Exp,Opt_))/(2*h);
  end
  % HStrain contribution, added in quadrature
  wH2 = 0;
  if isfield(Sys,'HStrain')
    [~,~,WH] = resfreqs_matrix(Sys,Exp,Opt_);
    wH2 = WH.^2;
  end
  Wfd = sqrt(sum((J*Cov).*J,2) + wH2);
  ok(k) = areequal(W,Wfd,1e-4,'rel');
end

end

%-------------------------------------------------------------------------------
function [f,idx] = parseref(str,Sys)
tok = regexp(str,'^(\w+)(\(.*\))?$','tokens','once');
f = tok{1};
if numel(tok)<2 || isempty(tok{2})
  idx = 1;
else
  sub = str2num(tok{2}(2:end-1)); %#ok<ST2NM>
  if isscalar(sub), idx = sub; else, idx = sub2ind(size(Sys.(f)),sub(1),sub(2)); end
end
end
