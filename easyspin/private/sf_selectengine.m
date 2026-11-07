% sf_selectengine  Choose the simulation engine of saffron
%
%   [method,reason,seq,reasons] = sf_selectengine(Sys,Exp,Opt)
%
% Checks which features of the request the subspace engine supports and
% returns the engine to use:
%   method   'subspace' or 'fullspace'
%   reason   why the full-space engine is used ('' for subspace)
%   seq      parsed sequence for the subspace engine (empty for fullspace)
%   reasons  cell array of the features the subspace engine does not
%            support (empty if the full-space engine was requested)
%
% Opt.Method = 'auto' (default) uses the subspace engine whenever it
% supports every feature, and the full-space engine otherwise. An explicit
% 'subspace' errors with the list of unsupported features.
%
% Capability table: the subspace engine supports one electron spin with
% S = 1/2 (S > 1/2 only when requested explicitly, as a transition-selective
% approximation), ideal and rectangular microwave pulses, single rectangular
% RF pulses with an RF frequency sweep, linear delay increments in the binning
% schemes, other pulse parameters point by point, coherence pathway
% selection, phase cycles, scalar T1/T2 and single-point S+ detection at
% the end of the sequence. Everything else needs the full-space engine (see
% sf_parsesequence for the sequence features). Conversely, inhomogeneous
% broadening (Sys.lwAnisoG) and ideal pulses with an excitation window
% (Opt.IdealPulses) are only supported by the subspace engine.

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
  [seq,reasons] = sf_parsesequence(Exp,Opt);
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
    reasons{end+1} = 'initial states (Sys.initState)';
  end
  if isfield(S_,'eqState') && ~isempty(S_.eqState)
    reasons{end+1} = 'equilibrium states (Sys.eqState)';
  end
  if (isfield(S_,'T1') && numel(S_.T1)>1) || (isfield(S_,'T2') && numel(S_.T2)>1)
    reasons{end+1} = 'transition-selective relaxation times (matrices in Sys.T1, Sys.T2)';
  end
  if isfield(S_,'nn') && any(S_.nn(:)~=0)
    reasons{end+1} = 'nuclear-nuclear couplings (Sys.nn)';
  end
  if ~isfield(S_,'Nucs') || isempty(S_.Nucs)
    reasons{end+1} = 'spin systems without nuclei';
  end
end
if highSpin && ~strcmp(requested,'subspace')
  reasons{end+1} = 'electron spin S > 1/2 (use Opt.Method = ''subspace'' for the transition-selective approximation)';
end
if isfield(Opt,'Relaxation') && numel(Opt.Relaxation)>1
  reasons{end+1} = 'event-specific relaxation (Opt.Relaxation)';
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
