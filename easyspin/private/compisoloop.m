% Provides a centralized loop over components and isotopologues that all
% simulation functions (pepper/garlic/chili/salt/saffron) need.
%
% Inputs:
%   simfc               simulation function (@pepper, @garlic, etc
%   Sys                 spin system(s)
%   Exp                 experimental parameters
%   Opt                 simulation options
%   autoRange           whether the sweep range should be determined automatically
%   thirdOutput         whether the third output argument has been requested 
%   separateComponents  whether to return components/isotopologues separately
%
% Outputs:
%   x      abscissa
%   spec   simulated spectrum or time trace
%   info   structure with additional information

function [x,spec,info] = compisoloop(simfcn,Sys,Exp,Opt,autoRange,thirdOutput,separateComponentSpectra)

includeInverseDomain = isequal(simfcn,@saffron);

info = struct;

if ~iscell(Sys), Sys = {Sys}; end

nComponents = numel(Sys);
logmsg(1,'  number of component spin systems: %d',nComponents);

% Determine isotopologues for each component
for c = 1:nComponents
  SysList{c} = isotopologues(Sys{c},Opt.IsoCutoff);  %#ok
  nIsotopologues(c) = numel(SysList{c});  %#ok
  logmsg(1,'    component %d: %d isotopologues',c,nIsotopologues(c));
end
nTotalComponents = sum(nIsotopologues);

if autoRange
  if nTotalComponents>1
    if Exp.FrequencySweep
      str = 'Exp.mwRange or Exp.mwCenterSweep';
    else
      str = 'Exp.Range or Exp.CenterSweep';
    end
    error('For multiple components, EasySpin cannot automatically determine a sweep range.\n Please specify sweep range manually using %s.',str);
  end
end

if separateComponentSpectra
  spec = [];
  if includeInverseDomain
    data_invdomain = [];
  end
else
  spec = 0;
  if includeInverseDomain
    data_invdomain = 0;
  end
end

% Initialize output structure
if thirdOutput
  info = struct;
  info.Component = [];
  info.Isotopologue = [];
  if any(nIsotopologues>1)
    info.Nucs = {};
  end
  info.Transitions = {};
  idx = 0;
end

% Loop over all components and isotopologues
for iComponent = 1:nComponents

  for iIsotopologue = 1:nIsotopologues(iComponent)

    % Simulate single-isotopologue spectrum
    Sys_ = SysList{iComponent}(iIsotopologue);
    Sys_.singleiso = true;
    [x,spec_,info_] = simfcn(Sys_,Exp,Opt);
    fdProvided = isfield(info_,'fd');

    % Accumulate or append spectra
    if separateComponentSpectra
      %if isvector(spec_), catdim = 1; else, catdim = ndims(spec_)+1; end
      catdim = 1;
      spec = cat(catdim,spec,spec_*Sys_.weight);
      if includeInverseDomain && fdProvided
        data_invdomain = cat(catdim,data_invdomain,info_.fd*Sys_.weight);
      end
    else
      spec = spec + spec_*Sys_.weight;
      if includeInverseDomain && fdProvided
        data_invdomain = data_invdomain + info_.fd*Sys_.weight;
      end
    end

    if thirdOutput
      idx = idx + 1;
      info.Component = [info.Component iComponent];
      info.Isotopologue = [info.Isotopologue iIsotopologue];
      if any(nIsotopologues>1)
        if isfield(Sys_,'Nucs')
          info.Nucs{idx} = Sys_.Nucs;
        else
          info.Nucs{idx} = '';
        end
      end
      if isfield(info_,'Transitions')
        info.Transitions{idx} = info_.Transitions;
        info.nTransitions(idx) = size(info_.Transitions,1);
      end
      if isfield(info_,'nSites')
        info.nSites(idx) = info_.nSites;
      end
      if isfield(info_,'nOrientations')
        info.nOrientations(idx) = info_.nOrientations;
      end
      if isfield(info_,'resfields')
        info.resfields{idx} = info_.resfields;
      end
    end

  end
end

if thirdOutput
  if numel(info.Transitions)==1
    info.Transitions = info.Transitions{1};
  end
  if isempty(info.Transitions)
    info = rmfield(info,'Transitions');
  end
  if isfield(info,'resfields') && numel(info.resfields)==1
    info.resfields = info.resfields{1};
  end
  if includeInverseDomain % for saffron
    if isfield(info_,'SinglePointDetection') % SimulationMode = 'thyme'
      info.SinglePointDetection = info_.SinglePointDetection;
    else  % SimulationMode = 'fast'
      if ~isempty(info_.td)
        info_.td = spec;
      end
      if fdProvided
        info_.fd = data_invdomain;
      end
      info_fields = fieldnames(info_);
      for i = 1:numel(info_fields)
        info.(info_fields{i}) = info_.(info_fields{i});
      end
    end
  end
end

end
