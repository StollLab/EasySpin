% orca2easyspin_maintxt   Read EPR properties from main ORCA output file

function [Sys,info] = orca2easyspin_maintxt(mainfile)

% File import and checks
%--------------------------------------------------------------------------
% Read entire file into cell array
fh = fopen(mainfile);
allLines = textscan(fh,'%s','Whitespace','','Delimiter','\n');
fclose(fh);
L = allLines{1};

% Remove empty lines and lines containing a single (non-printable) character
rmv = cellfun(@(x)length(strtrim(x))<=1,L);
L(rmv) = [];
nLines = numel(L);

% Assert that this is an ORCA output file
isOrcaOutputFile = any(contains(L{2},'* O   R   C   A *'));
if ~isOrcaOutputFile
  error('This is not an ORCA output file.');
end

% Determine ORCA version
versionLine = contains(L,'Program Version');
OrcaVersion = regexp(L{versionLine},'\d+\.\d+\.\d+','match','once');
info.OrcaVersion = OrcaVersion;


% Extract contents of input file
%--------------------------------------------------------------------------
k = 1;
while L{k}(1)~='|', k = k+1; end
startInput = k;
while L{k}(1)=='|', k = k+1; end
endInput = k - 2;
inputFile = L(startInput:endInput);
inputFile = regexprep(inputFile,'^\|\s*\d+>\s+','');
info.InputFile = char(inputFile);


% Determine whether the file contains multiple structures
%--------------------------------------------------------------------------
RunType(1).Header = '* Single Point Calculation *';
RunType(1).Step = '';
RunType(2).Header = '* Multiple XYZ Scan Calculation *'; % < v4
RunType(2).Step = 'MULTIPLE XYZ STEP'; % < v4
RunType(3).Header = '* Parameter Scan Calculation *';
RunType(3).Step = 'TRAJECTORY STEP';
RunType(4).Header = '*    Relaxed Surface Scan    *';
RunType(4).Step = 'RELAXED SURFACE SCAN STEP';

% Look for overall header
type = 0;
while k<=nLines && type==0
  for iType = 1:4
    if contains(L{k},RunType(iType).Header)
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
  startIdx = find(contains(L,RunType(type).Step));
  nStructures = numel(startIdx);
else
  startIdx = k;
  nStructures = 1;
end


% Loop over all structures and read properties
%--------------------------------------------------------------------------
data = struct;
for iStructure = 1:nStructures

  if iStructure<nStructures
    krange = startIdx(iStructure):startIdx(iStructure+1)-1;
  else
    krange = startIdx(iStructure):nLines;
  end

  % Atom info, and Cartesian coordinates
  %------------------------------------------------------------------------
  k = findheader('CARTESIAN COORDINATES (ANGSTROEM)',L,krange);
  k = k + 2;
  iAtom = 0;
  clear xyz Element NucId
  while L{k}(1)~='-'
    iAtom = iAtom + 1;
    xyz(iAtom,:) = sscanf(L{k},'%*s %f %f %f');  %#ok<AGROW>
    Element{iAtom} = sscanf(L{k},'%s',1);  %#ok<AGROW>
    NucId(iAtom) = elementsymbol2no(Element{iAtom});  %#ok<AGROW>
    k = k+1;
  end
  nAtoms = size(xyz,1);
  data(iStructure).nAtoms = nAtoms;
  data(iStructure).NucId = NucId;
  data(iStructure).Element = Element;
  data(iStructure).xyz = xyz;

  % Total charge
  %------------------------------------------------------------------------
  for k = krange
    if regexp(L{k},'^\s*Total Charge','match','once')
      break
    end
  end
  Charge = str2double(regexp(L{k},'\d+$','match','once'));
  data(iStructure).Charge = Charge;

  % Spin multiplicity
  %------------------------------------------------------------------------
  for k = krange
    if regexp(L{k},'^\s*Multiplicity','match','once')
      break
    end
  end
  Multiplicity = str2double(regexp(L{k},'\d+$','match','once'));
  data(iStructure).Multiplicity = Multiplicity;
  data(iStructure).S = (Multiplicity-1)/2;
  
  % Mulliken atomic charges and spin populations
  %------------------------------------------------------------------------
  MullikenTitle{1} = 'MULLIKEN ATOMIC CHARGES AND SPIN DENSITIES';  % <2.7
  MullikenTitle{2} = 'MULLIKEN ATOMIC CHARGES AND SPIN POPULATIONS'; % >=2.7
  MullikenSection = false;
  for k = krange
    if strcmp(L{k},MullikenTitle{1}) || strcmp(L{k},MullikenTitle{2})
      MullikenSection = true;
      break
    end
  end
  if MullikenSection
    k = k+2;
    Mulliken = zeros(nAtoms,2);
    for iAtom = 1:nAtoms
      Mulliken(iAtom,:) = sscanf(L{k+iAtom-1}(9:end),'%f %f').';
    end
  else
    Mulliken = zeros(0,2);
    fprintf('No Mulliken analysis found for structure %d of %d.\n',iStructure,nStructures);
  end
  data(iStructure).MullikenCharge = Mulliken(:,1);
  data(iStructure).MullikenSpin = Mulliken(:,2);

  % g matrix
  %------------------------------------------------------------------------
  k = findheader('ELECTRONIC G-MATRIX',L,krange);
  if ~isempty(k)
    k = k+3;
    % read raw asymmetric g matrix
    g_raw = readmatrix(L(k:k+2));
    g_sym = (g_raw.'*g_raw)^(1/2);
    g_sym = (g_sym+g_sym.')/2; % symmetrize numerically

    [V,g] = eig(g_sym);
    gvals = diag(g).';
    if det(V)<0
      V(:,1) = -V(:,1);
    end
    gFrame = eulang(V);
  else
    g_raw = [];
    g_sym = [];
    gvals = [];
    gFrame = [];
  end
  data(iStructure).graw = g_raw;
  data(iStructure).g = g_sym;
  data(iStructure).gvals = gvals;
  data(iStructure).gFrame = gFrame;

  % Zero-field splitting tensor
  %------------------------------------------------------------------------
  k = findheader('ZERO-FIELD-SPLITTING TENSOR',L,krange);
  if ~isempty(k)
    % read raw D matrix (cm^-1) and diagonalize
    D_raw = readmatrix(L(k+3:k+5));
    recalcVecs = true;
    if recalcVecs
      [D_vecs,D_vals] = eig(D_raw);
      D_vals = diag(D_vals).';
      if det(D_vecs)<0
        D_vecs(:,1) = -D_vecs(:,1);
      end
    else
      D_vals = sscanf(L{k+8},'%f %f %f').';
      D_vecs = readmatrix(L(k+9:k+11));
    end
    D_raw = D_raw*100*clight/1e6;  % cm^-1 -> MHz
    D_vals = D_vals*100*clight/1e6;  % cm^-1 -> MHz
    DFrame = eulang(D_vecs);
  else
    D_raw = [];
    D_vals = [];
    DFrame = [];
  end
  data(iStructure).Draw = D_raw;
  data(iStructure).Dvals = D_vals;
  data(iStructure).DFrame = DFrame;

  % Hyperfine and electric field gradient
  %------------------------------------------------------------------------
  Araw = cell(1,nAtoms);
  efg = cell(1,nAtoms);
  Q = cell(1,nAtoms);
  QFrame = cell(1,nAtoms);
  A = cell(1,nAtoms);
  AFrame = cell(1,nAtoms);

  k = findheader('ELECTRIC AND MAGNETIC HYPERFINE STRUCTURE',L,krange);
  if ~isempty(k)
    iAtom = 0;
    Araw = cell(1,nAtoms);
    efg = cell(1,nAtoms);
    while k<=krange(end)
      if regexp(L{k},'^\s*Nucleus\s*')
        iAtom = sscanf(L{k}(9:end),'%d',1)+1;
        [~,qrefEl] = referenceisotope(Element{iAtom});
      elseif regexp(L{k},'^\s*Raw HFC matrix')
        if strncmp(L{k+1}(2:4),'---',3)
          idx = k+2;
        else
          idx = k+1;
        end
        % Read raw HFC matrix
        Araw_ = readmatrix(L(idx:idx+2));
        Araw{iAtom} = Araw_;
        % Move to line starting with A(Tot) (version-dependent)
        idx = idx+3;
        while L{idx}(1:7)~=" A(Tot)"
          idx = idx+1;
        end
        % Read principal values and eigenvectors
        A_vals = sscanf(L{idx}(13:end),'%f %f %f').';
        R = readmatrix(L(idx+2:idx+4),5);
        A{iAtom} = A_vals;
        AFrame{iAtom} = eulang(R);
        k = k+5;

      elseif regexp(L{k},'^\s*Raw EFG matrix\s*')
        if strncmp(L{k+1}(2:4),'---',3)
          idx = k+2;
        else
          idx = k+1;
        end
        efg_ = readmatrix(L(idx:idx+2));
        efg{iAtom} = efg_; % atomic unit (Eh/e/a0^2)

        if ~isempty(qrefEl)  % element has I>=1 isotopes
          % Get EFG tensor and diagonalize
          efg_ = efg{iAtom}*hartree/echarge/bohrrad^2; % atomic unit -> SI unit
          [R,eq] = eig(efg_);
          eq = diag(eq).';

          % Sort by eigenvalue magnitude (this is the common convention)
          [~,idx] = sort(abs(eq));
          eq = eq(idx);
          R = R(:,idx);
          if det(R)<0
            R(:,1) = -R(:,1);
          end

          % Calculate quadrupole parameters and quadrupole tensor
          Qmom = qrefEl.qm*barn;
          I = qrefEl.I;
          e2qQh = echarge*Qmom*eq(3)/planck/1e6;
          K = e2qQh/(4*I)/(2*I-1);
          eta = (eq(1)-eq(2))/eq(3);
          Q{iAtom} = K*[-(1-eta),-(1+eta),2];
          QFrame{iAtom} = eulang(R);
        end
        k = k+3;
      end
      k = k+1;
    end
  end

  data(iStructure).hfc = Araw;
  data(iStructure).efg = efg;
  data(iStructure).Q = Q;
  data(iStructure).QFrame = QFrame;
  data(iStructure).A = A;
  data(iStructure).AFrame = AFrame;

end % for iStructure = 1:nStructures

% Copy relevant data to spin system structure
%--------------------------------------------------------------------------
for iStructure = nStructures:-1:1
  d = data(iStructure);
  
  % Coordinates
  if ~isempty(d.xyz)
    Sys(iStructure).xyz= d.xyz;
  end

  % Spin multiplicity
  Sys(iStructure).S = d.S; 

  % g tensor
  if ~isempty(d.g)
    Sys(iStructure).g = d.gvals;
  end
  if ~isempty(d.gFrame)
    Sys(iStructure).gFrame = d.gFrame;
  end

  % D tensor
  if ~isempty(d.Dvals)
    Sys(iStructure).D = d.Dvals;
  end
  if ~isempty(d.DFrame)
    Sys(iStructure).DFrame = d.DFrame;
  end

  % Compile nuclear data (isotopes, hyperfine coupling, quuadropole coupling)
  idx = 0;
  for iAtom = 1:nAtoms
    if ~isempty(d.A{iAtom}) || ~isempty(d.Q{iAtom})
      idx = idx + 1;
      Sys(iStructure).Nucs{idx} = d.Element{iAtom};
      if ~isempty(d.A{iAtom})
        Sys(iStructure).A(idx,1:3) = d.A{iAtom};
        Sys(iStructure).AFrame(idx,1:3) = d.AFrame{iAtom};
      end
      if ~isempty(d.Q{iAtom})
        Sys(iStructure).Q(idx,1:3) = d.Q{iAtom};
        Sys(iStructure).QFrame(idx,1:3) = d.QFrame{iAtom};
      else
        if isfield(Sys(iStructure),'Q')
          Sys(iStructure).Q(idx,1:3) = 0;
          Sys(iStructure).QFrame(idx,1:3) = 0;
        end
      end
    end
  end

  % Store all other data in spin system structure
  Sys(iStructure).data = data;
end

end


function M = readmatrix(L,startidx)
if nargin<2, startidx = 1; end
M(1,:) = sscanf(L{1}(startidx:end),'%f %f %f').';
M(2,:) = sscanf(L{2}(startidx:end),'%f %f %f').';
M(3,:) = sscanf(L{3}(startidx:end),'%f %f %f').';
end

function k = findheader(header,L,krange)
header_found = false;
for k = krange
  if strcmp(L{k},header)
    header_found = true;
    break
  end
end
if ~header_found
  k = [];
end
end
