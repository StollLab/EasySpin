% isotopologues  generates list of isotopologues
%
%    isotopologues(NucList)
%    isotopologues(NucList,Abundances)
%    out = ...
%    
%    Computes all possible isotopologues of the
%    given list of nuclei or elements, including
%    their natural abundances.
%
%    NucList: list of isotopes or elements
%      e.g. 'Cu', '63Cu,N,N', 'H,H,C,C'
%
%    out: structure containing the isotope lists
%      for the various isotopologues and their abundances
%
%    If out is not requested, the results are display.

function varargout = isotopologues(NucList,Abundances,nEquiv,IsoThreshold)

%------------------------------------------------------------------
if (nargin<4), IsoThreshold = 0.001; end
if (nargin<3), nEquiv = []; end
Debug = 0;
AbundanceSort = 1;
ConsolidateNonmagneticIsotopes = 1;
%------------------------------------------------------------------

if nargin==0, help(mfilename); return; end

if (nargin==1)
  Abundances = [];
end

if isempty(NucList)
  out.Nucs = {''};
  out.Abund = 1;
  out.nIso = 1;
  varargout = {out};
  return
end

if isstruct(NucList)
  if isfield(NucList,'Nucs')
    NucList = NucList.Nucs;
  else
    error('Please give list of nuclei as first argument (e.g. ''Cu,Cl,H,H'').');
  end
end

global IsotopeList % initialized by nucdata()
if isempty(IsotopeList)
  % if not initialized yet, call nucdata
  dummy = nucdata('1H');
end

if Debug, fprintf('Input string:  %s\n', NucList); end
NucList = nucstringparse(NucList,'m');
if Debug, disp('Initial parsing:'); NucList, end

if ~isempty(Abundances)
  if ~iscell(Abundances), Abundances = {Abundances}; end
  if numel(Abundances)<numel(NucList)
    error('Please specify all abundances!');
  end
end

if isempty(nEquiv), nEquiv = ones(1,numel(NucList)); end
if numel(nEquiv)~=numel(NucList)
  error('Number of isotopes and list of equivalent nuclei do not match.');
end

nNucs = numel(NucList);
Nucs = cell(1,nNucs);
Abund = cell(1,nNucs);

for iNuc = 1:numel(NucList)
  N = NucList{iNuc};
  
  if isempty(N)
    error('Missing nucleus in list, maybe a comma too many?');
  end
  
  if Debug, fprintf('===== List entry %d:  %s\n',iNuc,N); end
 
  if N(1)=='('
    
    if isempty(Abundances)
      error('For custom isotope mixtures, abundances must be given in Sys.Abund!');
    end
    
    if Debug, disp('Multiple isotopes'); end
    idx = strfind(N,')');
    ElementSymbol = N(idx+1:end);
    MassNumbers = sscanf(N(2:idx-1),'%d,');
    for k = 1:numel(MassNumbers)
      Nucs{iNuc}{k} = sprintf('%d%s',MassNumbers(k),ElementSymbol);
    end
    ab = Abundances{iNuc};
    Abund{iNuc} = ab/sum(ab);
    
  else
    
    MassNumber = sscanf(N,'%d');
    if isempty(MassNumber)
      % no mass number given ->
      % expand into natural abundance mixture
      if Debug, disp('Natural abundance expansion'); end
      
      iidx = find(strcmp(N,IsotopeList.Element));
      % compile list of available isotopes
      if ~isempty(iidx)
        MixData = [IsotopeList.Nucleons(iidx) IsotopeList.Abundances(iidx)];
      else
        error(sprintf('Could not find element ''%s''',N));
      end
      MixData(MixData(:,2)==0,:) = [];
      for k=1:size(MixData,1)
        Nucs{iNuc}{k} = sprintf('%d%s',MixData(k,1),N); 
      end
      Abund{iNuc} = MixData(:,2).'/sum(MixData(:,2));
     
    else
      
      if Debug, disp('Single isotope'); end
      Nucs{iNuc} = {N};
      Abund(iNuc) = {1};
      
    end
  end
end

if Debug
  disp('== Complete list ==');
  for iNuc=1:numel(Nucs)
    Nucs{iNuc}
    Abund{iNuc}
  end
end

nNucs = numel(Nucs);

if ~isempty(nEquiv)
  if numel(nEquiv)~=nNucs
    error('nEquiv array and number of nuclei are not consistent.');
  end
  for k=1:nNucs
    if numel(Nucs{k})>1 && nEquiv(k)>1
      error('Cannot compute isotope patterns for equivalent nuclei.');
    end
  end
end


% ConsolidateNonmagneticIsotopes list if wanted
if ConsolidateNonmagneticIsotopes
  if Debug, disp('Consolidating non-magnetic isotopes.'); end
  for k=1:numel(Nucs)
    idx = find(nucspin(Nucs{k})==0);
    if ~isempty(idx)
      Abund{k}(idx(1)) = sum(Abund{k}(idx));
      Nucs{k}(idx(2:end)) = [];
      Abund{k}(idx(2:end)) = [];
    end
  end
  if Debug
    disp('== Complete consolidated list ==');
    for k=1:numel(Nucs)
      Nucs{k}
      Abund{k}
    end
  end
end


% Automatic hyperfine and quadrupole reference isotope
%-------------------------------------------------------
for iNuc = nNucs:-1:1

  abundance_ = Abund{iNuc};
  [I,gn,qm] = nucdata(Nucs{iNuc});
  if isempty(I)
    error('Isotope with unknown spin encountered.');
  end
  
  abundance_(I<1/2) = -1;
  [mx,idx] = max(abundance_);
  gnref(iNuc) = gn(idx);
  gnscale_{iNuc} = gn/gnref(iNuc);
  
  abundance_(I<1) = -1;
  [mx,idx] = max(abundance_);
  if (mx>0)
    Qref(iNuc) = idx;
    qmref(iNuc) = qm(idx);
    qmscale_{iNuc} = qm/qmref(iNuc);
  else
    Qref(iNuc) = 0;
    qmref(iNuc) = NaN;
    qmscale_{iNuc} = ones(1,numel(qm));
  end
end

if Debug
  gnref
  Qref
  qmref
end

maxAbundance = 1;
for iNuc=1:nNucs
  nIsotopes(iNuc) = numel(Nucs{iNuc});
  maxAbundance = maxAbundance*max(Abund{iNuc});
end
minAbundance = IsoThreshold*maxAbundance;

nTotalIsotopologues = prod(nIsotopes);
if Debug
  fprintf('%d isotopologues total\n',nTotalIsotopologues);
end

% Abundance tree traversal:
% Compute list of isotopologues with their abundances
%--------------------------------------------------------------------------
iNuc = 1;
iIso = zeros(1,nNucs);
abundance = zeros(1,nNucs);
nIsotopologues = 0;
isoIdx = zeros(100,nNucs);
while (iNuc>=1)
  if (iNuc>nNucs) % if vertex, process and backtrack
    %if Debug
    %  fprintf('%7.6f  ',abundance(nNucs));
    %  for k=1:nNucs, fprintf('%d ',iIso(k)); end
    %  fprintf('\n');
    %end
    nIsotopologues = nIsotopologues + 1;
    if mod(nIsotopologues,100)==1
      isoIdx(nIsotopologues+100,:)=0;
    end
    isoIdx(nIsotopologues,:) = iIso;
    IsoAbund(nIsotopologues) = abundance(nNucs);
    iNuc = iNuc - 1;
  end
  if iIso(iNuc)<nIsotopes(iNuc) % next isotope
    iIso(iNuc) = iIso(iNuc) + 1;
    if iNuc>1, abund_ = abundance(iNuc-1); else abund_ = 1; end
    abund_ = abund_*Abund{iNuc}(iIso(iNuc));
    if abund_>=minAbundance % go to next nucleus
      abundance(iNuc) = abund_;
      iNuc = iNuc + 1;
    else % reset isotope index and go back to previous nucleus
      iIso(iNuc) = 0;
      iNuc = iNuc - 1;
    end
  else % reset isotope index and go back to previous nucleus
    iIso(iNuc) = 0;
    iNuc = iNuc - 1;
  end
end
isoIdx(nIsotopologues+1:end,:) = [];

% Compile isotope strings for isotopologues
%-------------------------------------------------------
% make string with comma-separated symbols
for k = 1:nIsotopologues
  str = [];
  for iNuc = 1:nNucs
    idx_ = isoIdx(k,iNuc);
    nu = Nucs{iNuc}{idx_};
    %ss{iNuc} = nu;
    str = [str ',' nu];
    Ascale{k}(iNuc) = gnscale_{iNuc}(idx_);
    Qscale{k}(iNuc) = qmscale_{iNuc}(idx_);
  end
  IsoNucsS{k} = str(2:end);
end


% Sort isotopologues 
%-------------------------------------------------------
if AbundanceSort
  [IsoAbund,idx] = sort(-IsoAbund); IsoAbund = -IsoAbund;
  IsoNucsS = IsoNucsS(idx);
  Ascale = Ascale(idx);
  Qscale = Qscale(idx);
end

if Debug
  for k=1:nIsotopologues
    fprintf('%0.7f  %s\n',IsoAbund(k),IsoNucsS{k});
  end
end

out.Nucs = IsoNucsS;
out.Abund = IsoAbund;
out.Ascale = Ascale;
out.Qscale = Qscale;
out.nIso = length(IsoNucsS);

% Print out
%-------------------------------------------------------
if (nargout==0)
  fprintf('Abs.Abund   Rel.Abund   Composition\n');
  for k=1:nIsotopologues
    fprintf(' %0.6f    %0.6f    %s\n',IsoAbund(k),IsoAbund(k)/max(IsoAbund),IsoNucsS{k});
  end
  fprintf('%d out of %d isotopologues with abundance above threshold\n',...
    nIsotopologues,nTotalIsotopologues);
else
  varargout = {out};
end
