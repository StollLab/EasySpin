% Generates list of components (specied and isotopologues)
% Sys: spin system structure, or cell array of spin system structures
% AbundanceCutoff: (optional) isotopologue abundance cutoff

function [SysList,weight] = expandcomponents(Sys,AbundanceCutoff)

if (nargin<2), AbundanceCutoff = 0; end

if isstruct(Sys), Sys = {Sys}; end

idx = 1;
for iSpecies = 1:numel(Sys)
  Sys_ = Sys{iSpecies};
  if ~isfield(Sys_,'weight')
    speciesWeight = 1;
  else
    speciesWeight = Sys_.weight;
  end
  if ~isfield(Sys_,'Nucs'), Sys_.Nucs = ''; end
  if ~isfield(Sys_,'Abund'), Sys_.Abund = []; end
  if ~isfield(Sys_,'n'), Sys_.n = []; end
  isoList = isotopologues(Sys_.Nucs,Sys_.Abund,Sys_.n,AbundanceCutoff);
  
  nIsotopologues = isoList.nIso;
  for iIsotopologue = 1:nIsotopologues
    Sys_.Nucs = isoList.Nucs{iIsotopologue};
    if ~isempty(Sys_.Nucs)
      Sys_.Ascale = isoList.Ascale{iIsotopologue};
      Sys_.Qscale = isoList.Qscale{iIsotopologue};
    end
    Sys_.singleiso = 1;
    
    SysList{idx} = Sys_;
    weight(idx) = speciesWeight*isoList.Abund(iIsotopologue);
    idx = idx+1;
  end % iIso
end % iComponent
