% isotopologues   Generate list of isotopologues
%
%    isotopologues(NucList)
%    isotopologues(NucList,nEquiv)
%    isotopologues(NucList,nEquiv,Abundances)
%    isotopologues(NucList,nEquiv,Abundances,relThreshold)
%    isotopologues(Sys)
%    isotopologues(Sys,relThreshold)
%    out = ...
%
%    Computes all possible isotopologues of the given list of nuclei or
%    elements, including their abundances.
%
%    Input:
%      NucList       list of isotopes or elements
%                       e.g. 'Cu', '63Cu,N,N', 'H,H,C,C'
%      nEquiv        number of equivalent nuclei (default: 1 for each
%                       nucleus in NucList)
%                       e.g. [1 4] for 'Cu,N'
%      Abundances    cell array of nuclear abundances
%      relThreshold  isotopologue abundance threshold, relative to
%                       abundance of most abundant isotopologue
%                       (default 0.001)
%
%    out                 structure array containing a list of all isotopologues
%       out(k).Nucs      string with list of isotopes
%       out(k).Abund     overall absolute abundance
%       out(k).n         number of equivalent nuclei
%
%    If out is not requested, the list of isotopologues is displayed.

function varargout = isotopologues(varargin)

if nargin==0, help(mfilename); return; end

SysInput = isstruct(varargin{1});

% Internal settings
%------------------------------------------------------------------
%#ok<*AGROW>
%#ok<*UNRCH>
defaultAbundanceThreshold = 1e-4;
consolidateNonmagneticIsotopes = true;
Debug = false;


% Parameters
%------------------------------------------------------------------
if SysInput
  Sys = varargin{1};
  relAbundanceThreshold = defaultAbundanceThreshold;
  if nargin==2
    relAbundanceThreshold = varargin{2};
  elseif nargin>2
    error('At most two inputs possible when a spin system is given.');
  end
  
  NucList = '';
  if isfield(Sys,'Nucs')
    NucList = Sys.Nucs;
  end
  nEquiv = [];
  if isfield(Sys,'n')
    nEquiv = Sys.n;
  end
  Abundances = [];
  if isfield(Sys,'Abund')
    Abundances = Sys.Abund;
  end
  
else
  Sys = [];
  NucList = varargin{1};
  nEquiv = [];
  Abundances = [];
  relAbundanceThreshold = defaultAbundanceThreshold;
  if nargin>=2, nEquiv = varargin{2}; end
  if nargin>=3, Abundances = varargin{3}; end
  if nargin==4, relAbundanceThreshold = varargin{4}; end
  if nargin>4
    error('At most four inputs are possible.');
  end
  
end

if ~ischar(NucList)
  if SysInput
    error('List of nuclei in Sys.Nucs must be a string (e.g. ''Cu,Cl,H,H'').');
  else
    error('List of nuclei must be a string (e.g. ''Cu,Cl,H,H'').');
  end
end


% Special case: no nuclei
%===============================================================================
if isempty(NucList)
  if nargout==0
    fprintf('Abs.Abund   Rel.Abund   Composition\n');
    fprintf('1.000       1.000       (no nuclei)\n');
  else
    if SysInput
      if ~isfield(Sys,'weight'), Sys.weight = 1; end
      isotopologue(1) = Sys;
    else
      isotopologue(1).Nucs = '';
      isotopologue(1).weight = 1;
      isotopologue(1).n = [];
    end
    varargout = {isotopologue};
  end
  return
end

global IsotopeList % initialized by nucdata()
if isempty(IsotopeList)
  % if not initialized yet, call nucdata()
  dummy = nucdata('1H');  %#ok<NASGU>
end

customMixtures = any(NucList=='(');
NucList = nucstring2list(NucList,'m');
nNucs = numel(NucList);

if customMixtures
  if isempty(Abundances)
    error('For custom isotope mixtures, abundances must be given!');
  end
end

if ~iscell(Abundances)
  if nNucs==1
    Abundances = {Abundances};
  end
  if isempty(Abundances)
    Abundances = cell(1,nNucs);
  end
else
  if numel(Abundances)~=nNucs
    error('Abundances must be a cell array with one entry (array) per nucleus.');
  end
end
if numel(Abundances)<nNucs
  error('%d nuclei, but only %d abundances given. Please specify all abundances!',...
    nNucs,numel(Abundances));
end

if isempty(nEquiv)
  nEquiv = ones(1,nNucs);
end
if numel(nEquiv)~=nNucs
  if SysInput
    error('%d nuclei given in Sys.Nucs, but %d numbers in Sys.n.',nNucs,numel(nEquiv));
  else
    error('%d nuclei given, but %d numbers in nEquiv.',nNucs,numel(nEquiv));
  end
end

Nucs = cell(1,nNucs);
Abund = cell(1,nNucs);

if Debug
  fprintf('Input:\n');
  fprintf('  NucList: ');
  for iNuc = 1:nNucs-1
    fprintf('%s,',NucList{iNuc});
  end
  fprintf('%s\n',NucList{end});
  fprintf('  nEquiv: %s\n',num2str(nEquiv));
end

% Determine format of hyperfine and quadrupole fields
if SysInput
  if isfield(Sys,'S')
    nElectrons = numel(Sys.S);
  else
    nElectrons = 1;
  end
  if isfield(Sys,'A')
    if nElectrons==1
      Aisotropic = numel(Sys.A)==nNucs;
    else
      Aisotropic = all(size(Sys.A)==[nNucs nElectrons]);
    end
    if Aisotropic
      Sys.A = reshape(Sys.A,nNucs,nElectrons);
    end
    Aaxial = all(size(Sys.A)==[nNucs 2*nElectrons]);
    Arhombic = all(size(Sys.A)==[nNucs 3*nElectrons]);
    Afull = all(size(Sys.A)==[3*nNucs 3*nElectrons]);
    Aexchange = size(Sys.A,1)==nNucs; % for compatiblity with chem. exchange program
  end
  if isfield(Sys,'A_')
    if nElectrons==1
      A_isotropic = numel(Sys.A_)==nNucs;
    else
      A_isotropic = all(size(Sys.A_)==[nNucs nElectrons]);
    end
    A_axial = all(size(Sys.A_)==[nNucs 2*nElectrons]);
    A_rhombic = all(size(Sys.A_)==[nNucs 3*nElectrons]);
    A_exchange = size(Sys.A_,1)==nNucs; % for compatiblity with chem. exchange program
  end
  if isfield(Sys,'Q')
    Qaxial = numel(Sys.Q)==nNucs;
    if Qaxial
      Sys.Q = Sys.Q(:);
    end
    Qrhombic = all(size(Sys.Q)==[nNucs 2]);
    Qpvalues = all(size(Sys.Q)==[nNucs 3]);
    Qfull = all(size(Sys.Q)==[3*nNucs 3]);
  end
end


% Process nuclear isotope string list
%===============================================================================
for iNuc = 1:nNucs
  
  % Parse nuclear list entry
  %-----------------------------------------------------------------------------
  % Input: NucList{iNuc}
  % Output: Element{iNuc}, isoAbundances
  N = NucList{iNuc};
  if isempty(N)
    error('Missing nucleus in list.');
  end
  if N(1)=='('
    % Mass numbers given -> custom-abundance isotope mixture
    
    idx = strfind(N,')');
    Element{iNuc} = N(idx+1:end);
    MassNumbers = sscanf(N(2:idx-1),'%d,');
    isoAbundances = Abundances{iNuc};
    
  elseif isempty(sscanf(N,'%d'))
    % No mass numbers given -> natural-abundance isotope mixture
    
    % Locate element and compile list of natural isotopes
    Element{iNuc} = N;
    iidx = find(strcmp(Element{iNuc},IsotopeList.Element));
    if ~isempty(iidx)
      MassNumbers = IsotopeList.Nucleons(iidx);
      isoAbundances = IsotopeList.Abundances(iidx).'/100;
    else
      error('Could not find element ''%s''',Element{iNuc});
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
    Element{iNuc} = N(floor(log10(MassNumbers)+1)+1:end);
    isoAbundances = 1;
    
  end
  Groups(iNuc).Element = Element{iNuc};
  Groups(iNuc).n = nEquiv(iNuc);
  
  % Get names and nuclear data for all isotopes
  %-----------------------------------------------------------------------------
  nIsotopes = numel(isoAbundances);
  for iIso = 1:nIsotopes
    Groups(iNuc).Isotopes{iIso} = sprintf('%d%s',MassNumbers(iIso),Element{iNuc});
  end
  [Groups(iNuc).I,gn,qm] = nucdata(Groups(iNuc).Isotopes);
  
  % Combine nonmagnetic isotopes if wanted
  %-----------------------------------------------------------------------------
  if consolidateNonmagneticIsotopes
    idx = find(Groups(iNuc).I==0);
    if ~isempty(idx)
      isoAbundances(idx(1)) = sum(isoAbundances(idx));
      rmv = idx(2:end);
      isoAbundances(rmv) = [];
      Groups(iNuc).Isotopes(rmv) = [];
      Groups(iNuc).I(rmv) = [];
      if SysInput
        gn(rmv) = [];
        qm(rmv) = [];
      end
    end
  end
  Groups(iNuc).Abundances = isoAbundances;
  
  % Calculate A and Q tensors for all isotopes
  %-----------------------------------------------------------------------------
  if SysInput
    
    % Determine hyperfine and quadrupole reference isotopes
    %------------------------------------------------------------
    abundances_ = Groups(iNuc).Abundances;
    % (a) find most abundant isotope with I>=1/2 for gn reference
    abundances_(Groups(iNuc).I<1/2) = -1;
    [mx,idx] = max(abundances_);
    if mx>0
      gnref = gn(idx);
    else
      gnref = [];
    end    
    % (b) find most abundant isotope with I>=1 for Q reference
    abundances_(Groups(iNuc).I<1) = -1;
    [mx,idx] = max(abundances_);
    if mx>0
      qmref = qm(idx);
    else
      qmref = [];
    end
    
    % Calculate A and Q tensors
    %------------------------------------------------------------
    % A - isotropic; axial; rhombic; full
    if isempty(gnref) || gnref==0, gnref = 1; end
    if isfield(Sys,'A')
      if Aisotropic || Aaxial || Arhombic || Aexchange
        for k = 1:numel(gn)
          Groups(iNuc).A{k} = Sys.A(iNuc,:)/gnref*gn(k);
        end
      elseif Afull
        for k = 1:numel(gn)
          Groups(iNuc).A{k} = Sys.A((iNuc-1)*3+(1:3),:)/gnref*gn(k);
        end
      else
        error('Sys.A size (%dx%d) is inconsistent with number of spins (%d electrons, %d nuclei).',...
          size(Sys.A,1),size(Sys.A,2),nElectrons,nNucs);
      end
    end
    
    % A_ - 1, 2, or 3 elements per nuclear spin
    if isfield(Sys,'A_')
      if A_isotropic
        for k = 1:numel(gn)
          Groups(iNuc).A_{k} = Sys.A_(iNuc)/gnref*gn(k);
        end
      elseif A_axial || A_rhombic
        for k = 1:numel(gn)
          Groups(iNuc).A_{k} = Sys.A_(iNuc,:)/gnref*gn(k);
        end
      else
        error('Sys.A_ size (%dx%d) is inconsistent with number of spins (%d electrons, %d nuclei).',...
          size(Sys.A_,1),size(Sys.A_,2),nElectrons,nNucs);
      end
    end
    
    % Q - Q alone; Q and eta; three principal values; full?
    if isempty(qmref) || qmref==0, qmref = 1; end
    if isfield(Sys,'Q')
      if Qaxial
        for k = 1:numel(qm)
          Groups(iNuc).Q{k} = Sys.Q(iNuc)/qmref*qm(k);
        end
      elseif Qrhombic || Qpvalues
        for k = 1:numel(gn)
          Groups(iNuc).Q{k} = Sys.Q(iNuc,:)/qmref*qm(k);
        end
      elseif Qfull
        for k = 1:numel(qm)
          Groups(iNuc).Q{k} = Sys.Q((iNuc-1)*3+(1:3),:)/qmref*qm(k);
        end
      else
        error('Sys.Q size (%dx%d) is inconsistent with number of nuclear spins (%d).',...
          size(Sys.Q,1),size(Sys.Q,2),nNucs);
      end
    end
    
    % AFrame and QFrame
    if isfield(Sys,'AFrame')
      if all(size(Sys.AFrame)==[nNucs,nElectrons*3])
        Groups(iNuc).AFrame = Sys.AFrame(iNuc,:);
      else
        error('Sys.AFrame has the wrong size.');
      end
    end
    if isfield(Sys,'QFrame')
      Groups(iNuc).QFrame = Sys.QFrame(iNuc,:);
    end
    
  end % if SysInput
  
end % for iNuc = 1:nNucs

if Debug
  disp('List of groups:')
  for k = 1:numel(Groups)
    fprintf('  Group %d\n',k);
    disp(Groups(k))
  end
end


% Generate all isotopologues within equivalent class
%===============================================================================
for iNuc = 1:nNucs
  nIsotopes = numel(Groups(iNuc).Isotopes);
  
  % Generate all isotopologues within equivalent class
  [kvec,multiplicity] = multisetlist(nEquiv(iNuc),nIsotopes);
  if any(isinf(multiplicity))
    error('Multiplicity overflow. - Cannot handle nEquiv=%d for nucleus #%d.',...
      nEquiv(iNuc),iNuc);
  end
  
  % Calculate symbol list, abundances, and #equivent nuclei
  for im = 1:numel(multiplicity)
    Nucs_ = {};
    abu_ = 1;
    equ_ = [];
    for iIso = 1:nIsotopes
      %if kvec(im,iIso)==0, continue; end
      abu_ = abu_*Groups(iNuc).Abundances(iIso)^kvec(im,iIso);
      Nucs_ = [Nucs_ Groups(iNuc).Isotopes(iIso)];
      equ_ = [equ_ kvec(im,iIso)];
    end
    if abu_==0
      error('Abundance underflow. - Cannot handle nEquiv=%d for nucleus #%d.',...
        nEquiv(iNuc),iNuc);
    end
    Nucs{iNuc}{im} = Nucs_;
    Abund{iNuc}(im) = abu_*multiplicity(im);
    nEquivList{iNuc}{im} = equ_;
  end
  if SysInput
    Groups(iNuc).kvec = kvec;
  end
end
maxAbundance = prod(cellfun(@max,Abund));
absAbundanceThreshold = relAbundanceThreshold*maxAbundance;
nIsotopes = cellfun(@numel,Nucs);


if Debug
  fprintf('Group isotopologues:\n');
  for iNuc = 1:numel(Nucs)
    fprintf('  Group %2d: %s\n',iNuc,Groups(iNuc).Element);
    for im = 1:numel(Nucs{iNuc})
      fprintf('    ');
      fprintf('%f  ',Abund{iNuc}(im));
      fprintf('%3d ',nEquivList{iNuc}{im});
      fprintf('\n');
    end
    fprintf('    sum of abundances: %f\n',sum(Abund{iNuc}));
  end
end


% Abundance tree traversal: calculate isotopologue abundances and isotope indices
%===============================================================================
[IsoListIdx,IsoListAbund] = abundancetreetraversal(nIsotopes,Abund,absAbundanceThreshold);
nIsotopologues = numel(IsoListAbund);


if Debug
  fprintf('List of isotopologues (abundances and indices)\n')
  for k = 1:nIsotopologues
    fprintf('  %f  ',IsoListAbund(k));
    for q = 1:size(IsoListIdx(k,:),2)
      fprintf('%2d ',IsoListIdx(k,q));
    end
    fprintf('\n');
  end
  fprintf('  %f   threshold\n',absAbundanceThreshold);
  fprintf('  sum of abundances: %f\n',sum(IsoListAbund));
end


% Compile isotopologue list
%===============================================================================
% make Nucs string with comma-separated isotope symbols
if SysInput
  for k = nIsotopologues:-1:1
    isotopologue(k) = Sys;
  end
end

for k = 1:nIsotopologues
  Nucs_ = {};
  n = [];
  A = [];
  A_ = [];
  Q = [];
  AFrame = [];
  QFrame = [];
  for iNuc = 1:nNucs
    idx = IsoListIdx(k,iNuc);
    nz = find(nEquivList{iNuc}{idx}~=0);
    gr = Groups(iNuc);
    for iIso = nz
      if gr.I(iIso)==0, continue; end
      Nucs_ = [Nucs_ Nucs{iNuc}{idx}(iIso)];
      n = [n nEquivList{iNuc}{idx}(iIso)];
      if SysInput
        if isfield(Sys,'A')
          if Aisotropic
            A = [A; gr.A{iIso}];
          elseif Aaxial || Arhombic || Aexchange || Afull
            A = [A; gr.A{iIso}];
          end
        end
        if isfield(Sys,'A_')
          if A_isotropic
            A_ = [A_ gr.A_{iIso}];
          elseif A_axial || A_rhombic
            A_ = [A_; gr.A_{iIso}];
          end
        end
        if isfield(Sys,'Q')
          if Qaxial || Qrhombic || Qpvalues || Qfull
            Q = [Q; gr.Q{iIso}];
          end
        end
        if isfield(Sys,'AFrame')
          AFrame = [AFrame; gr.AFrame];
        end
        if isfield(Sys,'QFrame')
          QFrame = [QFrame; gr.QFrame];
        end
      end
    end
  end
  
  isotopologue(k).Nucs = nuclist2string(Nucs_);
  isotopologue(k).n = n;
  
  if isfield(isotopologue(k),'weight') && ~isempty(isotopologue(k).weight)
    isotopologue(k).weight = isotopologue(k).weight*IsoListAbund(k);
  else
    isotopologue(k).weight = IsoListAbund(k);
  end
  
  if ~isempty(A)
    isotopologue(k).A = A;
  end
  if ~isempty(A_)
    isotopologue(k).A_ = A_;
  end
  if ~isempty(Q)
    isotopologue(k).Q = Q;
  end
  if ~isempty(AFrame)
    isotopologue(k).AFrame = AFrame;
  end
  if ~isempty(QFrame)
    isotopologue(k).QFrame = QFrame;
  end
  
  if isempty(isotopologue(k).Nucs)
    nucFields = {'A','A_','Q','AFrame','QFrame'};
    for iF = 1:numel(nucFields)
      if isfield(isotopologue(k),nucFields{iF})
        isotopologue(k).(nucFields{iF}) = [];
      end
    end
  end
  
end


if Debug
  fprintf('List of isotopologues\n');
  for k = 1:nIsotopologues
    fprintf(' isotopologue %d:\n',k);
    disp(isotopologue(k))
  end
  fprintf('\n');
end



% Output
%===============================================================================
if nargout==0
  
  % Determine longest isotope string
  strlen = 0;
  for k = 1:nIsotopologues
    strlen = max(strlen,length(isotopologue(k).Nucs));
  end
  
  % Print table of isotopologues
  fprintf('Sys.Abund   Sys.Nucs%sSys.n\n',repmat(' ',1,strlen-4));
  for k = 1:nIsotopologues
    isot = isotopologue(k);
    N = isot.Nucs;
    fprintf('  %0.6f   %s%s    ',...
      isot.weight,N,repmat(' ',1,strlen-length(N)));
    fprintf('%2d ',isot.n);
    fprintf('\n');
  end
  
  fprintf('%d of %d isotopologues with abundance above threshold (%0.2e)\n',...
    nIsotopologues,prod(nIsotopes),relAbundanceThreshold);
else
  varargout = {isotopologue};
end


%-------------------------------------------------------------------------------
%{
multisetlist   Generates list of isotopologue groups for a given number of
               of equivalent positions n and number of isotopes k.
 [kvec,multiplicity] = multisetlist(n,k)

 E.g. nIsotopes = k = 3 and nEquivPos = n = 4 yields kvec with 15 rows, one
 row for each group (multiset) of spectrally indistinguishable isotopologues.

 kvec =
     4     0     0
     3     1     0
     3     0     1
     2     2     0
     2     1     1
     2     0     2
     1     3     0
     1     2     1
     1     1     2
     1     0     3
     0     4     0
     0     3     1
     0     2     2
     0     1     3
     0     0     4

The meaning is the following. E.g., the second row states that this group of
isotopologues has the first isotope in 3 of the equivalent positions, the
second isotope in 1 position, and the third isotope in none.

multiplicity is a vector containing the multiplicities for each isotologue
group, i.e. the number of possible permutations among the positions that give
isotopologues that are indistinguishable. This number needs to be included in
the overall weight of a spectrum simulation.

For the example above,

multiplicity =
     1
     4
     4
     6
    12
     6
     4
    12
    12
     4
     1
     4
     6
     4
     1

The sum over all multiplicities equals k^n.
%}
%-------------------------------------------------------------------------------
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
  if j<k && v(j)>0
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


%-------------------------------------------------------------------------------
% Abundance tree traversal: Compute list of isotopologues with their abundances
%-------------------------------------------------------------------------------
function [IsoListIdx,IsoListAbund] = abundancetreetraversal(nIsotopes,Abund,absAbundanceThreshold)

if absAbundanceThreshold<0
  error('Abundance threshold cannot be negative.');
end

nNucs = numel(Abund);
iIso = zeros(1,nNucs);
abundance = zeros(nNucs,1);

IsoListIdx = [];
IsoListAbund = [];

nIsotopologues = 0;
iNuc = 1;
while iNuc>=1
  
  % if last nucleus done, add to list (if above threshold), and go back to
  % previous nucleus
  if iNuc>nNucs
    if abundance(nNucs)>absAbundanceThreshold
      nIsotopologues = nIsotopologues + 1;
      IsoListIdx(nIsotopologues,:) = iIso;
      IsoListAbund(nIsotopologues) = abundance(nNucs);
    end
    iNuc = iNuc - 1;
  else
    % if not last nucleus, update abundance and proceed to next group
    if iIso(iNuc)<nIsotopes(iNuc)
      iIso(iNuc) = iIso(iNuc) + 1;
      if iNuc>1
        abund_ = abundance(iNuc-1);
      else
        abund_ = 1;
      end
      abundance(iNuc) = abund_*Abund{iNuc}(iIso(iNuc));
      iNuc = iNuc + 1;
    else % reset isotope index and go back to previous nucleus
      iIso(iNuc) = 0;
      iNuc = iNuc - 1;
    end
  end
end
IsoListAbund = IsoListAbund.';

return