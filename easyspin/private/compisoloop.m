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
  info.Transitions = {};
  info.Components = [];
  info.Isotopologues = [];
  idx = 0;
end

% Loop over all components and isotopologues
for iComponent = 1:nComponents

  for iIsotopologue = 1:nIsotopologues(iComponent)

    % Simulate single-isotopologue spectrum
    Sys_ = SysList{iComponent}(iIsotopologue);
    Sys_.singleiso = true;
    [x,spec_,out_] = simfcn(Sys_,Exp,Opt);
    fdProvided = isfield(out_,'fd');
    size(spec_)

    % Accumulate or append spectra
    if separateComponentSpectra
      %if isvector(spec_), catdim = 1; else, catdim = ndims(spec_)+1; end
      catdim = 1;
      spec = cat(catdim,spec,spec_*Sys_.weight);
      if includeInverseDomain && fdProvided
        data_invdomain = cat(catdim,data_invdomain,out_.fd*Sys_.weight);
      end
    else
      spec = spec + spec_*Sys_.weight;
      if includeInverseDomain && fdProvided
        data_invdomain = data_invdomain + out_.fd*Sys_.weight;
      end
    end

    if thirdOutput
      idx = idx + 1;
      info.Components = [info.Components iComponent];
      info.Isotopologues = [info.Isotopologues iIsotopologue];
      if isfield(out_,'Transitions')
        info.Transitions{idx} = out_.Transitions;
        info.nTransitions(idx) = size(out_.Transitions,1);
      end
      if isfield(out_,'nSites')
        info.nSites(idx) = out_.nSites;
      end
      if isfield(out_,'nOrientations')
        info.nOrientations(idx) = out_.nOrientations;
      end
    end

  end
end
if thirdOutput
  if numel(info.Transitions)==1
    info.Transitions = info.Transitions{1};
  end
end

end
