% isotopologues   Generate list of isotopologues
%
%    isotopologues(NucList)
%    isotopologues(NucList,nEquiv)
%    isotopologues(NucList,nEquiv,Abundances)
%    isotopologues(NucList,nEquiv,Abundances,relThreshold)
%    out = ...
%    
%    Computes all possible isotopologues of the given list of nuclei or
%    elements, including their abundances.
%
%    NucList       list of isotopes or elements
%                     e.g. 'Cu', '63Cu,N,N', 'H,H,C,C'
%    nEquiv        number of equivalent nuclei (default: 1 for each
%                     nucleus in NucList)
%    Abundances    cell array of nuclear abundances
%    relThreshold  isotopologue abundance threshold, relative to
%                     abundance of most abundante isotopologue (default 0.001)
%
%    out                 structure array containing the isotopologue data
%       out(k).Nucs      string with list of isotopes
%       out(k).Abund     overall absolute abundance
%       out(k).n         number of equivalent nuclei
%       out(k).Ascale    scaling factor for hyperfine couplings
%       out(k).Qscale    scaling factor for quadrupole couplings
%
%    If out is not requested, the list of isotopologues is displayed.

function varargout = isotopologues(NucList,nEquiv,Abundances,relAbundanceThreshold)

if nargin==0, help(mfilename); return; end

Debug = false;
%#ok<*AGROW>
%#ok<*UNRCH> 
 
% Parameters
%------------------------------------------------------------------
if nargin<2, nEquiv = []; end
if nargin<3, Abundances = []; end
if nargin<4, relAbundanceThreshold = 1e-4; end
sortByAbundances = true;
consolidateNonmagneticIsotopes = true;

if iscell(NucList) && numel(NucList)==1
  NucList = NucList{1};
end


% Special case: no nuclei
%------------------------------------------------------------------
if isempty(NucList)
  isotopologue(1).Nucs = '';
  isotopologue(1).Abund = 1;
  isotopologue(1).n = 1;
  isotopologue(1).Ascale = 1;
  isotopologue(1).Qscale = 1;
  if nargout==0
    fprintf('Abs.Abund   Rel.Abund   Composition\n');
    fprintf('1.000       1.000       ''''\n');
  else
    varargout = {isotopologue};
  end
  return
end

global IsotopeList % initialized by nucdata()
if isempty(IsotopeList)
  % if not initialized yet, call nucdata
  dummy = nucdata('1H');  %#ok<NASGU>
end

if ~ischar(NucList)
  error('Please give list of nuclei as first argument (e.g. ''Cu,Cl,H,H'').');
end
customMixtures = any(NucList=='(');
NucList = nucstring2list(NucList,'m');

if customMixtures
  if isempty(Abundances)
    error('For custom isotope mixtures, abundances must be given!');
  end
end

if ~iscell(Abundances)
  if numel(NucList)==1
    Abundances = {Abundances};
  end
  if isempty(Abundances)
    Abundances = cell(1,numel(NucList));
  end
else
  if numel(Abundances)~=numel(NucList)
    error('Abundances must be a cell array with one entry (array) per nucleus.');
  end
end
if numel(Abundances)<numel(NucList)
  error('%d nuclei, but only %d abundances given. Please specify all abundances!',...
    numel(NucList),numel(Abundances));
end

if isempty(nEquiv)
  nEquiv = ones(1,numel(NucList));
end
if numel(nEquiv)~=numel(NucList)
  error('%d nuclei given, but only %d elements in nEqui.',numel(NucList),numel(nEquiv));
end

nNucs = numel(NucList);
Nucs = cell(1,nNucs);
Abund = cell(1,nNucs);

if Debug
  fprintf('Input:\n');
  fprintf('  NucList: ');
  for iNuc = 1:numel(NucList)-1
    fprintf('%s,',NucList{iNuc});
  end
  fprintf('%s\n',NucList{end});
  fprintf('  nEquiv: %s\n',num2str(nEquiv));
end

% Process nuclear isotope string list
%-------------------------------------------------------------------------------
for iNuc = 1:numel(NucList)
  
  % Parse nuclear string entry
  N = NucList{iNuc};
  if isempty(N)
    error('Missing nucleus in list.');
  end  
  if N(1)=='('
    % Mass numbers given -> custom-abundance isotope mixture
    
    idx = strfind(N,')');
    ElementSymbol{iNuc} = N(idx+1:end);
    MassNumbers = sscanf(N(2:idx-1),'%d,');
    isoAbundances = Abundances{iNuc};
        
  elseif isempty(sscanf(N,'%d'))
    % No mass numbers given -> natural-abundance isotope mixture
    
    % Locate element and compile list of natural isotopes
    ElementSymbol{iNuc} = N;
    iidx = find(strcmp(ElementSymbol{iNuc},IsotopeList.Element));
    if ~isempty(iidx)
      MassNumbers = IsotopeList.Nucleons(iidx);
      isoAbundances = IsotopeList.Abundances(iidx)/100;
    else
      error('Could not find element ''%s''',ElementSymbol{iNuc});
    end
    idx = isoAbundances==0;
    MassNumbers(idx) = [];
    isoAbundances(idx) = [];
    
    if ~isempty(Abundances{iNuc})
      error('Nucleus %d is a natural-abundance mixture; provide [] in Abundances{%d}.',iNuc,iNuc);
    end
  else
    % Explicit single isotope
    
    MassNumbers = sscanf(N,'%d');
    ElementSymbol{iNuc} = N(floor(log10(MassNumbers)+1)+1:end);
    isoAbundances = 1;
    
  end
  
  % Get names and nuclear data for all isotopes
  nIsotopes = numel(MassNumbers);
  for iIso = 1:nIsotopes
    IsoName{iNuc}{iIso} = sprintf('%d%s',MassNumbers(iIso),ElementSymbol{iNuc});
  end
  [I{iNuc},gn{iNuc},qm{iNuc}] = nucdata(IsoName{iNuc});
  
  % Combine nonmagnetic isotopes if wanted
  if consolidateNonmagneticIsotopes
    idx = find(I{iNuc}==0);
    if ~isempty(idx)
      isoAbundances(idx(1)) = sum(isoAbundances(idx));
      rmv = idx(2:end);
      isoAbundances(rmv) = [];
      MassNumbers(rmv) = [];
      IsoName{iNuc}(rmv) = [];
      I{iNuc}(rmv) = [];
      gn{iNuc}(rmv) = [];
      qm{iNuc}(rmv) = [];
    end
    nIsotopes = numel(MassNumbers);
  end
  IsotopeAbundances{iNuc} = isoAbundances;
  
  % Determine hyperfine and quadrupole reference isotopes
  abundances = IsotopeAbundances{iNuc};
  
  % (a) find most-abundant isotope with I>=1/2 and set its gn as reference
  abundances(I{iNuc}<1/2) = -1;
  [mx,idx] = max(abundances);
  if (mx>0)
    gnref(iNuc) = gn{iNuc}(idx);
    gnscale_{iNuc} = gn{iNuc}/gnref(iNuc);
  else
    gnref(iNuc) = NaN;
    gnscale_{iNuc} = ones(1,numel(gn{iNuc}));
  end
  
  % (b) find most-abundant isotope with I>=1 and set its Q as reference
  abundances(I{iNuc}<1) = -1;
  [mx,idx] = max(abundances);
  if (mx>0)
    qmref(iNuc) = qm{iNuc}(idx);
    qmscale_{iNuc} = qm{iNuc}/qmref(iNuc);
  else
    qmref(iNuc) = NaN;
    qmscale_{iNuc} = ones(1,numel(qm{iNuc}));
  end
  
  % Generate all isotopologues within equivalent class
  [kvec,multiplicity] = multisetlist(nEquiv(iNuc),nIsotopes);
  nSets = numel(multiplicity);
  if any(isinf(multiplicity))
    error('Multiplicity overflow. - Cannot handle nEquiv=%d for nucleus #%d.',...
      nEquiv(iNuc),iNuc);
  end
  
  % Calculate strings, abundances, and #equivent nuclei
  for im = 1:nSets
    str_ = '';
    abu_ = 1;
    equ_ = [];
    for iIso = 1:nIsotopes
      if kvec(im,iIso)==0, continue; end
      abu_ = abu_*isoAbundances(iIso)^kvec(im,iIso);
      str_ = [str_ ',' IsoName{iNuc}{iIso}];
      equ_ = [equ_ kvec(im,iIso)];
    end
    if (abu_==0)
      error('Abundance underflow. - Cannot handle nEquiv=%d for nucleus #%d.',...
        nEquiv(iNuc),iNuc);
    end
    if numel(str_)>1, str_ = str_(2:end); end
    Nucs{iNuc}{im} = str_;
    Abund{iNuc}(im) = abu_*multiplicity(im);
    nEquivList{iNuc}{im} = equ_;
  end
    kvecc{iNuc} = kvec;
end


if Debug
  
  fprintf('Parsed list:\n');
  for iNuc = 1:numel(Nucs)
    fprintf('  Nucleus group %2d: %s\n',iNuc,ElementSymbol{iNuc});
    for im = 1:numel(Nucs{iNuc})
      fprintf('    ');
      fprintf('%f  ',Abund{iNuc}(im));
      fprintf('%-15s  ',Nucs{iNuc}{im});
      fprintf('%3d ',nEquivList{iNuc}{im});
      fprintf('\n');
    end
    fprintf('    sum of abundances: %f\n',sum(Abund{iNuc}));
  end
  
end

maxAbundance = prod(cellfun(@max,Abund));
absAbundanceThreshold = relAbundanceThreshold*maxAbundance;
nIsotopes = cellfun(@numel,Nucs);
nTotalIsotopologues = prod(nIsotopes);


% Abundance tree traversal: Compute list of isotopologues with their abundances
%-------------------------------------------------------------------------------
iNuc = 1;
iIso = zeros(1,nNucs);
abundance = zeros(1,nNucs);
nIsotopologues = 0;
isoIdx = zeros(100,nNucs);
while (iNuc>=1)
  if (iNuc>nNucs) % if vertex, process and go to previous
    nIsotopologues = nIsotopologues + 1;
    if mod(nIsotopologues,100)==1 % allocate new memory every 100
      isoIdx(nIsotopologues+100,:) = 0;
    end
    isoIdx(nIsotopologues,:) = iIso;
    IsoAbund(nIsotopologues) = abundance(nNucs);
    iNuc = iNuc - 1;
  end
  if iIso(iNuc)<nIsotopes(iNuc) % next isotope
    iIso(iNuc) = iIso(iNuc) + 1;
    if iNuc>1, abund_ = abundance(iNuc-1); else, abund_ = 1; end
    abund_ = abund_*Abund{iNuc}(iIso(iNuc));
    if abund_>=absAbundanceThreshold % save abundance and go to next nucleus
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

if Debug
  fprintf('List of isotopologues (abundances and indices)\n')
  for k = 1:nIsotopologues
    fprintf('  %f  ',IsoAbund(k));
    for q = 1:size(isoIdx(k,:),2)
      fprintf('%2d ',isoIdx(k,q));
    end
    fprintf('\n');
  end
  fprintf('  %f   threshold\n',absAbundanceThreshold);
  fprintf('  sum of abundances: %f\n',sum(IsoAbund));
end

% Compile isotopologue list
%-------------------------------------------------------------------------------
% make string with comma-separated symbols
for k = nIsotopologues:-1:1
  str = [];
  equiv_ = [];
  Ascale_ = [];
  Qscale_ = [];
  for iNuc = 1:nNucs
    idx_ = isoIdx(k,iNuc);
    idx1 = find(kvecc{iNuc}(idx_,:));
    str = [str ',' Nucs{iNuc}{idx_}];
    equiv_ = [equiv_ nEquivList{iNuc}{idx_}];
    Ascale_ = [Ascale_ gnscale_{iNuc}(idx1)];
    Qscale_ = [Qscale_ qmscale_{iNuc}(idx1)];
  end
  isotopologue(k).Nucs = str(2:end); % exclude the initial comma
  isotopologue(k).Abund = IsoAbund(k);
  isotopologue(k).n = equiv_;
  isotopologue(k).Ascale = Ascale_;
  isotopologue(k).Qscale = Qscale_;
end

if Debug
end



% Sort isotopologues in descending order of abundance
%-------------------------------------------------------------------------------
if sortByAbundances && 0
  [dummy,idx] = sort([isotopologue.Abund],'descend'); %#ok<ASGLU>
  isotopologue = isotopologue(idx);
end


% Output
%-------------------------------------------------------------------------------
if (nargout==0)
  
  % Determine longest isotope string
  strlen = 0;
  for k = 1:nIsotopologues
    strlen = max(strlen,length(isotopologue(k).Nucs));
  end
  
  % Print table of isotopologues
  fprintf('abundance  Sys.Nucs%sSys.n\n',repmat(' ',1,strlen-4));  
  for k = 1:nIsotopologues
    isot = isotopologue(k);
    fprintf('  %0.6f   %s%s    ',...
      isot.Abund,isot.Nucs,repmat(' ',1,strlen-length(isot.Nucs)));
    fprintf('%2d ',isot.n);
    fprintf('\n');
  end
  
  fprintf('%d of %d isotopologues with abundance above threshold (%0.2e)\n',...
    nIsotopologues,nTotalIsotopologues,relAbundanceThreshold);
else
  varargout = {isotopologue};
end


%-------------------------------------------------------------------------------
% multisetlist   Generates list of isotopologue groups for a given number of
%                of equivalent positions n and number of isotopes k.
%
%  [kvec,multiplicity] = multisetlist(n,k)
%
%  E.g. nIsotopes = k = 3 and nEquivPos = n = 4 yields kvec with 15 rows, one
%  row for each group (multiset) of spectrally indistinguishable isotopologue.
%
%  kvec =
%      4     0     0
%      3     1     0
%      3     0     1
%      2     2     0
%      2     1     1
%      2     0     2
%      1     3     0
%      1     2     1
%      1     1     2
%      1     0     3
%      0     4     0
%      0     3     1
%      0     2     2
%      0     1     3
%      0     0     4
%
% The meaning is the following. E.g., the second row states that this group of
% isotopologues has the first isotope in 3 positions, the second isotope in 1
% position, and the third isotope in none.
%
% multiplicity a vector containing the multiplicities for each isotologue group,
% i.e. the number of possible permutations among the positions that give
% isotopologues that are indistinguishable. This number needs to be included in
% the overall weight of a spectrum simulation.
%
% For the example above,
%
% multiplicity =
%      1
%      4
%      4
%      6
%     12
%      6
%      4
%     12
%     12
%      4
%      1
%      4
%      6
%      4
%      1
%
% The sum over all multiplicities equals k^n.

function [kvec,multiplicity] = multisetlist(n,k)

% k = total number of different isotopes
% n = total number of equivalent positions

nMultiSets = nchoosek(n+k-1,k-1);
kvec = zeros(nMultiSets,k);
v = zeros(1,k);
v(1) = n;
idx = 1;
j = 1;
while j>0
  
  % store
  kvec(idx,:) = v;
  idx = idx + 1;
  
  % move right if possible
  if (j<k) && v(j)>0
    % push one right
    v(j) = v(j)-1;
    v(j+1) = v(j+1)+1;
    % step right
    j = j + 1;
  else
    % take value and step left as much as needed
    s = v(j);
    v(j) = 0;
    if s==n, break; end
    while v(j)==0
      j = j-1;
    end
    % push one right, and set number to the right
    v(j) = v(j)-1;
    v(j+1) = s+1;
    % step right
    j = j+1;
  end
  
end

multiplicity = zeros(nMultiSets,1);
for iSet = 1:nMultiSets
  % calculate multiplicity via multinomial coefficient
  k_ = [];
  for i = 1:k
    k_ = [k_ 1:kvec(iSet,i)];
  end
  multiplicity(iSet) = prod((1:n)./sort(k_));
end

return
