% orca2easyspin_maintxt   Read EPR properties from main ORCA output file

function [Sys,info] = orca2easyspin_maintxt(mainfile,HyperfineCutoff)

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
  CoordinateHeader = 'CARTESIAN COORDINATES (ANGSTROEM)';
  for k = krange
    if strcmp(L{k},CoordinateHeader)
      break
    end
  end
  k = k + 2;
  iAtom = 0;
  clear xyz Element NucId
  while L{k}(1)~='-'
    iAtom = iAtom + 1;
    xyz(iAtom,:) = sscanf(L{k},'%*s %f %f %f');
    Element{iAtom} = sscanf(L{k},'%s',1);
    NucId(iAtom) = elementsymbol2no(Element{iAtom});
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
  gMatrixHeader = 'ELECTRONIC G-MATRIX';
  gMatrixData = false;
  for k = krange
    if strcmp(L{k},gMatrixHeader)
      gMatrixData = true;
      break
    end
  end
  
  if gMatrixData
    k = k+3;
    % read raw asymmetric g matrix
    g_raw(1,:) = sscanf(L{k+0},'%f %f %f').';
    g_raw(2,:) = sscanf(L{k+1},'%f %f %f').';
    g_raw(3,:) = sscanf(L{k+2},'%f %f %f').';
    g_sym = (g_raw.'*g_raw)^(1/2);
    g_sym = (g_sym+g_sym.')/2; % symmetrize numerically

    [V,g] = eig(g_sym);
    if det(V)<0
      idx = [1 3 2];
      V = V(:,idx);
      g = g(idx,idx);
    end
    gvals = diag(g).';
    gFrame = eulang(V);
  else
    g_raw = [];
    g_sym = [];
    gvals = [];
    gFrame = [];
  end
  data(iStructure).g_raw = g_raw;
  data(iStructure).g = g_sym;
  data(iStructure).gvals = gvals;
  data(iStructure).gFrame = gFrame;

  % D matrix
  %------------------------------------------------------------------------
  DMatrixHeader = 'ZERO-FIELD-SPLITTING TENSOR';
  DMatrixData = false;
  for k = krange
    if strcmp(L{k},DMatrixHeader)
      DMatrixData = true;
      break
    end
  end

  if DMatrixData
    k = k+3;
    % read raw D matrix (cm^-1)
    D_raw(1,:) = sscanf(L{k+0},'%f %f %f').';
    D_raw(2,:) = sscanf(L{k+1},'%f %f %f').';
    D_raw(3,:) = sscanf(L{k+2},'%f %f %f').';

    [V,D] = eig(D_raw);
    D = diag(D).';
    if det(V)<0
      % make sure eigenvectors form right-handed coordinate frame
      idx = [1 3 2];
      V = V(:,idx);
      D = D(idx);
    end
    Dvals = D*100*clight/1e6;  % cm^-1 -> MHz
    DFrame = eulang(V);
  else
    D_raw = [];
    Dvals = [];
    DFrame = [];
  end

  data(iStructure).D_raw = D_raw;
  data(iStructure).Dvals = Dvals;
  data(iStructure).DFrame = DFrame;

  % Hyperfine and electric field gradient
  %-----------------------------------------------------------------
  hfc = cell(1,nAtoms);
  efg = cell(1,nAtoms);
  Q = cell(1,nAtoms);
  QFrame = cell(1,nAtoms);
  A = cell(1,nAtoms);
  AFrame = cell(1,nAtoms);

  nuclearHeader = 'ELECTRIC AND MAGNETIC HYPERFINE STRUCTURE';
  HyperfineEFGData = false;
  for k = krange
    if strcmp(L{k},nuclearHeader)
      HyperfineEFGData = true;
      break
    end
  end

  if HyperfineEFGData
    iAtom = 0;
    hfc = cell(1,nAtoms);
    efg = cell(1,nAtoms);
    while k<=krange(end)
      if regexp(L{k},'^\s*Nucleus\s*')
        iAtom = sscanf(L{k}(9:end),'%d',1)+1;
        [grefEl,qrefEl] = referenceisotope(Element{iAtom});
      elseif regexp(L{k},'^\s*Raw HFC matrix\s*')
        idx = k+1;
        if strncmp(L{k+1}(2:4),'---',3)
          idx = k+2;
        end
        hfc_(1,:) = sscanf(L{idx},'%f %f %f');
        hfc_(2,:) = sscanf(L{idx+1},'%f %f %f').';
        hfc_(3,:) = sscanf(L{idx+2},'%f %f %f').';
        hfc{iAtom} = hfc_;
        [R,A_] = eig(hfc_);
        A_ = diag(A_);
        if det(R)<0
          % make sure eigenvectors form right-handed coordinate frame
          idx = [1 3 2];
          R = R(:,idx);
          A_ = A_(idx);
        end
        A{iAtom} = A_;
        AFrame{iAtom} = eulang(R);
        k = k+3;
      elseif regexp(L{k},'^\s*Raw EFG matrix\s*')
        idx = k+1;
        if strncmp(L{k+1}(2:4),'---',3)
          idx = k+2;
        end
        efg_(1,:) = sscanf(L{idx},'%f %f %f');
        efg_(2,:) = sscanf(L{idx+1},'%f %f %f').';
        efg_(3,:) = sscanf(L{idx+2},'%f %f %f').';
        efg{iAtom} = efg_; % atomic unit

        efg_ = efg{iAtom}*hartree/echarge/bohrrad^2; % atomic unit -> SI unit
        [R,eq] = eig(efg_);
        if det(R)<0
          % make sure eigenvectors form right-handed coordinate frame
          idx = [1 3 2];
          R = R(:,idx);
          eq = eq(idx,idx);
        end
        eq = diag(eq);
        eta = (eq(1)-eq(2))/eq(3);
        if ~isempty(qrefEl)  % element has I>=1 isotopes
          Qmom = qrefEl.qm*barn;
          I = qrefEl.I;
          e2qQh = echarge*Qmom*eq(3)/planck/1e6;
          K = e2qQh/(4*I)/(2*I-1);
          Q{iAtom} = K*[-(1-eta),-(1+eta),2];
          QFrame{iAtom} = eulang(R);
        end
        k = k+3;
      end
      k = k+1;
    end
  end

  data(iStructure).hfc = hfc;
  data(iStructure).efg = efg;
  data(iStructure).Q = Q;
  data(iStructure).QFrame = QFrame;
  data(iStructure).A = A;
  data(iStructure).AFrame = AFrame;

end % for iStructure = 1:nStructures

% Copy relevant data to spin system structure
%--------------------------------------------------------------------------
Sys.S = data.S;
for iStructure = 1:nStructures
  d = data(iStructure);
  if ~isempty(d.g)
    Sys(iStructure).g = d.gvals;
  end
  if ~isempty(d.gFrame)
    Sys(iStructure).gFrame = d.gFrame;
  end
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
      end
    end
  end

  Sys(iStructure).data = data;
end

% Build spin system(s)
%--------------------------------------------------------------------------
%{
Sys = struct;
for iStructure = 1:nStructures
  if gMatrixData
    Sys(iStructure).g = data(iStructure).gvals;
    Sys(iStructure).gFrame = data(iStructure).gFrame;
  end
  if DMatrixData
    Sys(iStructure).D = data(iStructure).Dvals;
    Sys(iStructure).DFrame = data(iStructure).DFrame;
  end
  idx = 1;
  for k = 1:numel(hfc)
    A = data(iStructure).hfc{k};
    if isempty(A), continue; end
    [V,A] = eig(A); A = diag(A).';
    if det(V)<0, V = V(:,[1 3 2]); end
    Apa = eulang(V);
    efg = data(iStructure).efg{k};
    if ~isempty(efg)
      [V,efg] = eig(efg);
      efg = diag(efg).';
      if det(V)<0, V = V(:,[1 3 2]); end
      Qpa = eulang(V);
      eta = data(iStructure).eta(k);
      eeqQscaled = data(iStructure).eeqQscaled(k);
      P = eeqQscaled*[-(1-eta),-(1+eta),2];
      Sys(iStructure) = nucspinadd(Sys(iStructure),data(iStructure).Element{k},A,Apa,P,Qpa);
    else
      Sys(iStructure) = nucspinadd(Sys(iStructure),data(iStructure).Element{k},A,Apa);
    end
    idx = idx + 1;
  end
end

%}

end
