% sf_selectengine  Choose the simulation engine of saffron
%
%   [method,reason,seq,reasons] = sf_selectengine(Sys,Exp,Opt)
%
%   Checks which features of the request the subspace engine supports and
%   returns the engine to use. Opt.Method = 'auto' (default) uses the
%   subspace engine whenever it supports every feature, and the full-space
%   engine otherwise. An explicit 'subspace' errors with the list of
%   unsupported features.
%
%   Capability table: the subspace engine supports one electron spin with
%   S = 1/2 (S > 1/2 only when requested explicitly, as a transition-selective
%   approximation), ideal and rectangular microwave pulses, single rectangular
%   RF pulses with an RF frequency sweep, linear delay increments (also
%   moving pulses, binned or point by point, see Opt.Binning), other delay
%   and pulse parameter increments point by point, coherence pathway
%   selection, phase cycles, scalar T1/T2 (also only during some delays,
%   Opt.Relaxation), electron-only initial states (Sys.initState) on powders
%   without relaxation, and single-point S+ detection at the end of the
%   sequence or shifted from it (scalar Exp.DetWindow). Everything else needs the full-space engine (see
%   sf_parsesequence for the sequence features). Conversely, inhomogeneous
%   broadening (Sys.lwAnisoG) and ideal pulses with an excitation window
%   (Opt.IdealPulses) are only supported by the subspace engine.
%
%   Input:
%     Sys  spin system, or cell array of spin systems (components)
%     Exp  experiment, passed on to sf_parsesequence
%     Opt  options structure, with fields
%       Method       'auto' (default), 'subspace', 'fullspace'
%       Relaxation   relaxation on/off, one value or one per element of
%                    Exp.Sequence, passed on to sf_parsesequence
%       IdealPulses  true: treat finite pulses as ideal
%       crystalSample  true for crystals (from p_sampletype)
%
%   Output:
%     method   engine to use, 'subspace' or 'fullspace'
%     reason   why this engine is used, for the log; '' if the subspace
%              engine is selected automatically
%     seq      parsed sequence for the subspace engine (empty for fullspace)
%     reasons  cell array of the features the subspace engine does not
%              support; if the full-space engine was requested, the
%              sequence is not parsed and only spin system features are
%              listed

function [method,reason,seq,reasons] = sf_selectengine(Sys,Exp,Opt)

if ~isfield(Opt,'Method') || isempty(Opt.Method)
  Opt.Method = 'auto';
end
[~,err] = parseoption(Opt,'Method',{'auto','subspace','fullspace'});
error(err);
requested = char(Opt.Method);

seq = [];
reasons = {};
if ~strcmp(requested,'fullspace')
  idealPulses = isfield(Opt,'IdealPulses') && Opt.IdealPulses;
  relaxation = [];
  if isfield(Opt,'Relaxation'), relaxation = Opt.Relaxation; end
  [seq,reasons] = sf_parsesequence(Exp,idealPulses,relaxation);
end

% Spin system features
if iscell(Sys), SysList = Sys; else, SysList = {Sys}; end
highSpin = false;
for c = 1:numel(SysList)
  S_ = SysList{c};
  if isfield(S_,'S'), Svals = S_.S; else, Svals = 1/2; end
  if numel(Svals)>1
    reasons{end+1} = 'more than one electron spin'; %#ok<*AGROW>
  elseif Svals~=1/2
    highSpin = true;
  end
  if isfield(S_,'initState') && ~isempty(S_.initState)
    % Electron-only states (dimension prod(2S+1)) on powders, without
    % relaxation; states with a nuclear part, and crystals, need the
    % full-space engine (the shorthands such as 'singlet' are for two
    % electron spins, a reason of its own)
    state = S_.initState;
    if iscell(state), state = state{1}; end
    if isnumeric(state) && max(size(state))~=prod(2*Svals+1)
      reasons{end+1} = 'initial states with a nuclear part (Sys.initState)';
    end
    if isfield(Opt,'crystalSample') && Opt.crystalSample
      reasons{end+1} = 'initial states in crystals (Sys.initState)';
    end
    relaxationTimes = (isfield(S_,'T1') && any(S_.T1(:))) || (isfield(S_,'T2') && any(S_.T2(:)));
    relaxationOn = ~isfield(Opt,'Relaxation') || isempty(Opt.Relaxation) || any(Opt.Relaxation);
    if relaxationTimes && relaxationOn
      reasons{end+1} = 'initial states with relaxation (Sys.initState)';
    end
  end
  if isfield(S_,'eqState') && ~isempty(S_.eqState)
    reasons{end+1} = 'equilibrium states (Sys.eqState)';
  end
  if (isfield(S_,'T1') && numel(S_.T1)>1) || (isfield(S_,'T2') && numel(S_.T2)>1)
    reasons{end+1} = 'transition-selective relaxation times (matrices in Sys.T1, Sys.T2)';
  end
end
if highSpin && ~strcmp(requested,'subspace')
  reasons{end+1} = 'electron spin S > 1/2 (use Opt.Method = ''subspace'' for the transition-selective approximation)';
end
reasons = unique(reasons,'stable');

switch requested
  case 'fullspace'
    method = 'fullspace';
    reason = 'requested (Opt.Method)';
  case 'subspace'
    if ~isempty(reasons)
      error('Opt.Method = ''subspace'' cannot be used: the subspace engine does not support %s. Use ''fullspace'', or omit Opt.Method.',...
        strjoin(reasons,', '));
    end
    method = 'subspace';
    reason = 'requested (Opt.Method)';
  case 'auto'
    if isempty(reasons)
      method = 'subspace';
      reason = '';
    else
      method = 'fullspace';
      reason = ['the subspace engine does not support ' strjoin(reasons,', ')];
    end
end
if ~strcmp(method,'subspace'), seq = []; end

% Increments of delays before the first pulse have no effect on the initial
% state of the subspace engine (populations of the electron levels)
if strcmp(method,'subspace') && ~isempty(seq.prePulseDims)
  error('Exp.Dim%d: incrementing a delay before the first pulse has no effect.',seq.prePulseDims(1));
end

% Features only the subspace engine supports: inhomogeneous broadening of
% the EPR line (Sys.lwAnisoG, or its deprecated name Sys.HStrain)
if strcmp(method,'fullspace')
  for c = 1:numel(SysList)
    S_ = SysList{c};
    if (isfield(S_,'lwAnisoG') && any(S_.lwAnisoG(:))) || (isfield(S_,'HStrain') && any(S_.HStrain(:)))
      if strcmp(requested,'fullspace')
        error('The full-space engine (Opt.Method = ''fullspace'') does not support Sys.lwAnisoG.');
      else
        error(['Sys.lwAnisoG is only supported by the subspace engine, but this simulation needs ' ...
          'the full-space engine (%s).'],reason);
      end
    end
  end
end

% Ideal pulses with an excitation window (Opt.IdealPulses) are only
% supported by the subspace engine
if strcmp(method,'fullspace') && isfield(Opt,'IdealPulses') && Opt.IdealPulses
  if strcmp(requested,'fullspace')
    error('The full-space engine (Opt.Method = ''fullspace'') does not support Opt.IdealPulses = true.');
  else
    error(['Opt.IdealPulses = true is only supported by the subspace engine, but this simulation needs ' ...
      'the full-space engine (%s).'],reason);
  end
end

end
