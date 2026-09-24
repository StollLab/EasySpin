% orca2easyspin_maintxt   Read EPR properties from main ORCA output file

function [Sys,data] = orca2easyspin_mainout(mainfile)

% File import and checks
%--------------------------------------------------------------------------
% Read entire file into cell array
L = cellstr(readlines(mainfile));

% Remove empty lines & lines with a single character (to shorten file)
rmv = cellfun(@(x)length(strtrim(x))<=1,L);
L(rmv) = [];
nLines = numel(L);

% Assert that this is an ORCA output file
isOrcaOutputFile = nLines>=2 && contains(L{2},'* O   R   C   A *');
if ~isOrcaOutputFile
  error('This is not an ORCA output file.');
end

% Determine ORCA version
versionLine = find(contains(L,'Program Version'),1);
if isempty(versionLine)
  error('Could not find ORCA version information in the file.');
end
OrcaVersion = regexp(L{versionLine},'\d+\.\d+\.\d+','match','once');
data.OrcaVersion = OrcaVersion;

% Extract contents of input file
%--------------------------------------------------------------------------
k = 1;
while k<=nLines && L{k}(1)~='|', k = k+1; end
if k>nLines, error('Could not find start of input file echo.'); end
startInput = k;
while k<=nLines && L{k}(1)=='|', k = k+1; end
if k>nLines, error('Could not find end of input file echo.'); end
endInput = k - 2;
inputFile = L(startInput:endInput);
inputFile = regexprep(inputFile,'^\|\s*\d+>\s+','');
data.InputFile = char(inputFile);


% Determine whether the file contains multiple structures
%--------------------------------------------------------------------------
runType(1).Header = '* Single Point Calculation *';
runType(1).Step = '';
runType(2).Header = '* Multiple XYZ Scan Calculation *'; % < v4
runType(2).Step = 'MULTIPLE XYZ STEP'; % < v4
runType(3).Header = '* Parameter Scan Calculation *';
runType(3).Step = 'TRAJECTORY STEP';
runType(4).Header = '*    Relaxed Surface Scan    *';
runType(4).Step = 'RELAXED SURFACE SCAN STEP';

% Look for overall header
type = 0;
while k<=nLines && type==0
  for iType = 1:4
    if contains(L{k},runType(iType).Header)
      type = iType;
      break
    end
  end
  k = k + 1;
end
if k>nLines
  error('Could not determine type of run (single points, relaxed scan, etc.) from file.')
end

% Look for step titles
multipleStructures = type~=1;
if multipleStructures
  startIdx = find(contains(L,runType(type).Step));
  nStructures = numel(startIdx);
else
  startIdx = k;
  nStructures = 1;
end


% Loop over all structures and read properties
%--------------------------------------------------------------------------
data = repmat(data,1,nStructures);
for s = 1:nStructures

  if s<nStructures
    linerange = startIdx(s):startIdx(s+1)-1;
  else
    linerange = startIdx(s):nLines;
  end

  % Atom info and Cartesian coordinates
  [nAtoms,NucId,Element,xyz] = parsecoordinates(L,linerange);
  data(s).nAtoms = nAtoms;
  data(s).NucId = NucId;
  data(s).Element = Element;
  data(s).xyz = xyz;

  % Total charge
  charge = parsecharge(L,linerange);
  data(s).Charge = charge;

  % Spin multiplicity
  [Multiplicity,S] = parsespinmultiplicity(L,linerange);
  data(s).Multiplicity = Multiplicity;
  data(s).S = S;

  % Mulliken atomic charges and spin populations
  [MullikenChargePop,MullikenSpinPop] = parsemulliken(L,linerange,nAtoms);
  data(s).MullikenCharge = MullikenChargePop;
  data(s).MullikenSpin = MullikenSpinPop;

  % g matrix
  [graw,g,gvals,gFrame] = parsegmatrix(L,linerange);
  data(s).graw = graw;
  data(s).g = g;
  data(s).gvals = gvals;
  data(s).gFrame = gFrame;

  % Zero-field splitting tensor
  [Draw,Dvals,DFrame] = parsezfs(L,linerange);
  data(s).Draw = Draw;
  data(s).Dvals = Dvals;
  data(s).DFrame = DFrame;

  % Hyperfine and electric field gradient
  [Araw,Avals,AFrame,efg,Qvals,QFrame] = parsehyperfinequadrupole(L,linerange,Element,nAtoms);
  data(s).Araw = Araw;
  data(s).A = Avals;
  data(s).AFrame = AFrame;  
  data(s).efg = efg;
  data(s).Q = Qvals;
  data(s).QFrame = QFrame;

end  % for s = 1:nStructures

% Copy relevant data to spin system structure
%--------------------------------------------------------------------------
for s = nStructures:-1:1
  d = data(s);
  
  % Coordinates
  if ~isempty(d.xyz)
    Sys(s).xyz= d.xyz;
  end

  % Spin multiplicity, charge, elements
  Sys(s).S = d.S;
  Sys(s).charge = d.Charge;
  Sys(s).Elements = d.Element;

  % g tensor
  if ~isempty(d.g)
    Sys(s).g = d.gvals;
  end
  if ~isempty(d.gFrame)
    Sys(s).gFrame = d.gFrame;
  end

  % D tensor
  if ~isempty(d.Dvals)
    Sys(s).D = d.Dvals;
  end
  if ~isempty(d.DFrame)
    Sys(s).DFrame = d.DFrame;
  end

  % Compile nuclear data (isotopes, hyperfine coupling, quadrupole coupling)
  idx = 0;
  for iAtom = 1:d.nAtoms
    if ~isempty(d.A{iAtom}) || ~isempty(d.Q{iAtom})
      idx = idx + 1;
      Sys(s).Nucs{idx} = d.Element{iAtom};
      Sys(s).NucsIdx(idx) = iAtom;
      if ~isempty(d.A{iAtom})
        Sys(s).A(idx,1:3) = d.A{iAtom};
        Sys(s).AFrame(idx,1:3) = d.AFrame{iAtom};
      end
      if ~isempty(d.Q{iAtom})
        Sys(s).Q(idx,1:3) = d.Q{iAtom};
        Sys(s).QFrame(idx,1:3) = d.QFrame{iAtom};
      else
        if isfield(Sys(s),'Q')
          Sys(s).Q(idx,1:3) = 0;
          Sys(s).QFrame(idx,1:3) = 0;
        end
      end
    end
  end
  if idx>0
    Sys(s).Nucs = nuclist2string(Sys(s).Nucs);
  end

  % Store all other data in spin system structure
  Sys(s).data = data;
end

end

%-------------------------------------------------------------------------------
function M = parsematrix(L,startidx)
if nargin<2, startidx = 1; end
for k = 3:-1:1
  M(k,:) = sscanf(L{k}(startidx:end),'%f %f %f').';
end
end

%-------------------------------------------------------------------------------
function k = findheader(header,L,krange)
header_found = false;
for k = krange
  if strncmp(L{k},header,length(header))
    header_found = true;
    break
  end
end
if ~header_found
  k = [];
end
end

%-------------------------------------------------------------------------------
function [vals,angles] = diagonalizetensor(T,sortByMagnitude)
if nargin<2
  sortByMagnitude = false;
end
[R_T2M,T_] = eig(T);
vals = diag(T_).';
if sortByMagnitude
  [~,idx] = sort(abs(vals));
  vals = vals(idx);
  R_T2M = R_T2M(:,idx);
end
% Enforce right-handed frame
if det(R_T2M)<0
  R_T2M(:,1) = -R_T2M(:,1);
end
angles = eulang(R_T2M.');
end

%-------------------------------------------------------------------------------
function [nAtoms,NucId,Element,xyz] = parsecoordinates(L,krange)
k = findheader('CARTESIAN COORDINATES (ANGSTROEM)',L,krange);
if isempty(k)
  error('Cartesian coordinates not found.');
end
k = k + 2;
xyz = [];
Element = {};
NucId = [];
iAtom = 0;
while L{k}(1)~='-'
  iAtom = iAtom + 1;
  xyz(iAtom,:) = sscanf(L{k},'%*s %f %f %f');  %#ok<AGROW>
  Element{iAtom} = sscanf(L{k},'%s',1);  %#ok<AGROW>
  NucId(iAtom) = elementsymbol2no(Element{iAtom});  %#ok<AGROW>
  k = k+1;
end
nAtoms = size(xyz,1);
end

%-------------------------------------------------------------------------------
% In multi-structure files, the charge might be printed only once (e.g.
% for the first step of a parameter scan). If it is not found within
% krange, the entire file is searched.
function charge = parsecharge(L,krange)
k = findline(L,krange,'^\s*Total Charge');
if isempty(k)
  k = findline(L,1:numel(L),'^\s*Total Charge');
end
if isempty(k)
  error('Charge not found.');
end
charge = str2double(regexp(L{k},'-?\d+$','match','once'));
end

%-------------------------------------------------------------------------------
% Same fallback to the entire file as for the charge.
function [Multiplicity,S] = parsespinmultiplicity(L,krange)
k = findline(L,krange,'^\s*Multiplicity');
if isempty(k)
  k = findline(L,1:numel(L),'^\s*Multiplicity');
end
if isempty(k)
  error('Spin multiplicity not found.');
end
Multiplicity = str2double(regexp(L{k},'\d+$','match','once'));
S = (Multiplicity-1)/2;
end

%-------------------------------------------------------------------------------
% Index of first line in krange that matches the regular expression pattern
function k = findline(L,krange,pattern)
k = krange(find(~cellfun(@isempty,regexp(L(krange),pattern,'once')),1));
end

%-------------------------------------------------------------------------------
% The Mulliken analysis can be absent (e.g. with NoPop, or for later steps
% of a parameter scan). In that case, return empty arrays silently.
function [MullikenCharge,MullikenSpin] = parsemulliken(L,krange,nAtoms)
MullikenTitle{1} = 'MULLIKEN ATOMIC CHARGES AND SPIN DENSITIES';  % <2.7
MullikenTitle{2} = 'MULLIKEN ATOMIC CHARGES AND SPIN POPULATIONS'; % >=2.7
found = false;
for k = krange
  if strcmp(L{k},MullikenTitle{1}) || strcmp(L{k},MullikenTitle{2})
    found = true;
    break
  end
end
if found
  k = k+2;
  Mulliken = zeros(nAtoms,2);
  for iAtom = 1:nAtoms
    Mulliken(iAtom,:) = sscanf(L{k+iAtom-1}(9:end),'%f %f').';
  end
else
  Mulliken = zeros(0,2);
end
MullikenCharge = Mulliken(:,1);
MullikenSpin = Mulliken(:,2);
end

%-------------------------------------------------------------------------------
function [graw,gsym,gvals,gFrame] = parsegmatrix(L,searchrange)
k = findheader('ELECTRONIC G-MATRIX',L,searchrange);
if isempty(k)
  graw = [];
  gsym = [];
  gvals = [];
  gFrame = [];
  return
end

% Locate g matrix (number of lines down from header depends on ORCA version)
while ~contains(L{k},'The g-matrix')
  k = k+1;
end

% Read either the raw matrix & symmetrize & diagonalize it, or read the
% principal values and rotation matrix as calculated by ORCA. The latter is
% typically more accurate, since the former relies on the few sigfigs that
% are printed. But the difference might not matter.
readRawMatrix = false;
if readRawMatrix
  % Read raw (asymmetric) g matrix and symmetrize
  graw = parsematrix(L(k+(1:3)));
  gsym = sqrtm(graw*graw.');  % symmetrize
  gsym = (gsym+gsym.')/2;  % eliminate numerical errors
  % Diagonalize to get eigenvalues and orientation
  [gvals,gFrame] = diagonalizetensor(gsym);
else
  % Read symmetrized g principal values and rotation matrix
  while ~contains(L{k},'g(tot)')
    k = k+1;
  end
  gvals = sscanf(L{k}(12:end),'%f %f %f').';
  while ~contains(L{k},'Orientation:')
    k = k+1;
  end
  R_g2M = parsematrix(L(k+(1:3)),12);
  graw = [];
  R_M2g = R_g2M.';
  gFrame = eulang(R_M2g);
  gsym = R_g2M*diag(gvals)*R_g2M.'; 
end
end

%-------------------------------------------------------------------------------
function [Draw,Dvals,DFrame] = parsezfs(L,krange)
k = findheader('ZERO-FIELD-SPLITTING TENSOR',L,krange);
if isempty(k)
  Draw = [];
  Dvals = [];
  DFrame = [];
  return
end

readRawMatrix = true;
if readRawMatrix
  % Read raw D matrix (cm^-1) and convert to MHz
  Draw = parsematrix(L(k+3:k+5));
  Draw = Draw*100*clight/1e6;  % cm^-1 -> MHz
  % Diagonalize to get eigenvalues and orientation
  [Dvals,DFrame] = diagonalizetensor(Draw);
else
  % Read principal values and rotation matrix
  Dvals = sscanf(L{k+8},'%f %f %f').';
  Dvals = Dvals*100*clight/1e6;  % cm^-1 -> MHz
  R_D2M = parsematrix(L(k+8+(1:3)));
  DFrame = eulang(R_D2M);
  Draw = [];
end
end

function [Araw,A,AFrame,efg,Q,QFrame] = parsehyperfinequadrupole(L,krange,Element,nAtoms)

Araw = cell(1,nAtoms);
efg = cell(1,nAtoms);
Q = cell(1,nAtoms);
QFrame = cell(1,nAtoms);
A = cell(1,nAtoms);
AFrame = cell(1,nAtoms);

k = findheader('ELECTRIC AND MAGNETIC HYPERFINE STRUCTURE',L,krange);
if isempty(k)
  return
end

iAtom = 0;
qrefEl = [];
while k<=krange(end)
  if regexp(L{k},'ORCA EULER ANGLE PROGRAM'); break; end

  if regexp(L{k},'^\s*Nucleus\s+\d')
    iAtom = sscanf(L{k}(9:end),'%d',1)+1;
    [~,qrefEl] = referenceisotope(Element{iAtom});

  elseif regexp(L{k},'^\s*(Raw HFC matrix|Total HFC matrix)')
    if strncmp(L{k+1}(2:4),'---',3)
      idx = k+2;
    else
      idx = k+1;
    end
    % Read raw HFC matrix
    Araw{iAtom} = parsematrix(L(idx:idx+2));
    % Move to line starting with A(Tot) (version-dependent)
    idx = idx+3;
    while L{idx}(1:7)~=" A(Tot)"
      idx = idx+1;
    end
    % Read principal values and eigenvectors
    Avals = sscanf(L{idx}(13:end),'%f %f %f').';
    R_A2M = parsematrix(L(idx+2:idx+4),5);
    A{iAtom} = Avals;
    AFrame{iAtom} = eulang(R_A2M.');
    k = k+5;

  elseif regexp(L{k},'^\s*Raw EFG matrix\s*')
    if strncmp(L{k+1}(2:4),'---',3)
      idx = k+2;
    else
      idx = k+1;
    end
    efg{iAtom} = parsematrix(L(idx:idx+2)); % atomic unit (Eh/e/a0^2)

    if ~isempty(qrefEl)  % element has I>=1 isotopes
      % Get EFG tensor, diagonalize, sort by eigenvalue magnitude
      efg_SI = efg{iAtom}*hartree/echarge/bohrrad^2; % atomic unit -> SI unit
      [eq,QFrame{iAtom}] = diagonalizetensor(efg_SI,true);

      % Calculate quadrupole parameters and quadrupole tensor
      Qmom = qrefEl.qm*barn;
      I = qrefEl.I;
      e2qQh = echarge*Qmom*eq(3)/planck/1e6;
      K = e2qQh/(4*I)/(2*I-1);
      eta = (eq(1)-eq(2))/eq(3);
      Q{iAtom} = K*[-(1-eta),-(1+eta),2];
    end
    k = k+3;
  end
  k = k+1;
end

end
