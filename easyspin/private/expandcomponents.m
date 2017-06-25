% Generates list of components (species and isotopologues)
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
  isoList = isotopologues(Sys_.Nucs,Sys_.n,Sys_.Abund,AbundanceCutoff);
  
  nIsotopologues = numel(isoList);
  for k = 1:nIsotopologues
    Sys_.Nucs = isoList(k).Nucs;
    Sys_.n = isoList(k).n;
    if ~isempty(Sys_.Nucs)
      Sys_.Ascale = isoList(k).Ascale;
      Sys_.Qscale = isoList(k).Qscale;
    end
    Sys_.singleiso = 1;
    
    SysList{idx} = Sys_;
    weight(idx) = speciesWeight*isoList(k).Abund;
    idx = idx+1;
  end % iIso
end % iComponent
