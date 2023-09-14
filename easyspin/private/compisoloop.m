function [x,spec,out] = compisoloop(simfcn,Sys,Exp,Opt,autoRange,thirdOutput,separateComponentSpectra)

out = struct;

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
else
  spec = 0;
end

% Initialize output structure
if thirdOutput
  out = struct;
  out.Transitions = {};
  out.Components = [];
  out.Isotopologues = [];
  idx = 0;
end

% Loop over all components and isotopologues
for iComponent = 1:nComponents

  for iIsotopologue = 1:nIsotopologues(iComponent)

    % Simulate single-isotopologue spectrum
    Sys_ = SysList{iComponent}(iIsotopologue);
    Sys_.singleiso = true;
    [x,spec_,out_] = simfcn(Sys_,Exp,Opt);

    % Accumulate or append spectra
    if separateComponentSpectra
      spec = [spec; spec_*Sys_.weight];  %#ok
    else
      spec = spec + spec_*Sys_.weight;
    end

    if thirdOutput
      idx = idx + 1;
      out.Components = [out.Components iComponent];
      out.Isotopologues = [out.Isotopologues iIsotopologue];
      if isfield(out_,'Transitions')
        out.Transitions{idx} = out_.Transitions;
        out.nTransitions(idx) = size(out_.Transitions,1);
      end
      if isfield(out_,'nSites')
        out.nSites(idx) = out_.nSites;
      end
      if isfield(out_,'nOrientations')
        out.nOrientations(idx) = out_.nOrientations;
      end
    end

  end
end
if thirdOutput
  if numel(out.Transitions)==1
    out.Transitions = out.Transitions{1};
  end
end

end
