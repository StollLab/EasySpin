% orca2easyspin_maintxt   Read EPR properties from main ORCA output file

function [Sys,info] = orca2easyspin_maintxt(mainfile,HyperfineCutoff)

% File import and basic checks
%--------------------------------------------------------------------------
% Read entire file into cell array
fh = fopen(mainfile);
allLines = textscan(fh,'%s','whitespace','','delimiter','\n');
fclose(fh);
L = allLines{1};

% Remove empty line or those containing a single (non-printable) character
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


% Determine whether it's a multiple-structure file
%--------------------------------------------------------------------------
MultiXYZHeader{1} = '* Multiple XYZ Scan Calculation *'; % ORCA 3
multiStepTitle{1} = 'MULTIPLE XYZ STEP'; % ORCA 3
MultiXYZHeader{2} = '* Parameter Scan Calculation *';
multiStepTitle{2} = 'TRAJECTORY STEP';
MultiXYZHeader{3} = '*    Relaxed Surface Scan    *';
multiStepTitle{3} = 'RELAXED SURFACE SCAN STEP';
startIdx = k;
nStructures = 1;
multipleStructures = 0;
while k<=nLines && ~multipleStructures
  for q = 1:3
    if contains(L{k},MultiXYZHeader{q})
      multipleStructures = q;
      break
    end
  end
  k = k + 1;
end
if multipleStructures>0
  startIdx = find(contains(L,multiStepTitle{multipleStructures}));
  nStructures = numel(startIdx);
end


% Loop over all structures and read properties
%--------------------------------------------------------------------------
data = struct;
for iStructure = 1:nStructures

  kstart = startIdx(iStructure);
  if iStructure<nStructures
    kend = startIdx(iStructure+1)-1;
  else
    kend = nLines;
  end

  % Atom info, and Cartesian coordinates
  %------------------------------------------------------------------------
  CoordinateHeader = 'CARTESIAN COORDINATES (ANGSTROEM)';
  for k = kstart:kend
    if strcmp(L{k},CoordinateHeader)
      break;
    end
  end
  k = k + 2;
  iAtom = 1;
  clear xyz Element NucId
  while L{k}(1)~='-'
    xyz(iAtom,:) = sscanf(L{k},'%*s %f %f %f');
    Element{iAtom} = sscanf(L{k},'%s',1);
    NucId(iAtom) = elementsymbol2no(Element{iAtom});
    k = k+1;
    iAtom = iAtom + 1;
  end
  nAtoms = size(xyz,1);
  data(iStructure).nAtoms = nAtoms;
  data(iStructure).NucId = NucId;
  data(iStructure).Element = Element;
  data(iStructure).xyz = xyz;

  % Total charge
  %------------------------------------------------------------------------
  for k = kstart:kend
    if regexp(L{k},'^\s*Total Charge','match','once')
      break
    end
  end
  Charge = str2double(regexp(L{k},'\d+$','match','once'));
  data(iStructure).Charge = Charge;

  % Spin multiplicity
  %------------------------------------------------------------------------
  for k = kstart:kend
    if regexp(L{k},'^\s*Multiplicity','match','once')
      break
    end
  end
  Multiplicity = str2double(regexp(L{k},'\d+$','match','once'));
  data(iStructure).Multiplicity = Multiplicity;
  data(iStructure).S = (Multiplicity-1)/2;
  

  % Mulliken atomic charges and spin populations
  %------------------------------------------------------------------------
  % versions prior to 2.7
  MullikenHeader26 = 'MULLIKEN ATOMIC CHARGES AND SPIN DENSITIES';
  % version 2.7 and later
  MullikenHeader = 'MULLIKEN ATOMIC CHARGES AND SPIN POPULATIONS';
  MullikenSection = false;
  for k = kstart:kend
    if strcmp(L{k},MullikenHeader26) || strcmp(L{k},MullikenHeader)
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
  for k = kstart:kend
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
  for k = kstart:kend
    if strcmp(L{k},DMatrixHeader)
      DMatrixData = true;
      break
    end
  end

  if DMatrixData
    k = k+3;
    % read raw D matrix
    D_raw(1,:) = sscanf(L{k+0},'%f %f %f').';
    D_raw(2,:) = sscanf(L{k+1},'%f %f %f').';
    D_raw(3,:) = sscanf(L{k+2},'%f %f %f').';

    [V,D] = eig(D_raw);
    if det(V)<0
      % make sure eigenvectors form right-handed coordinate frame
      idx = [1 3 2];
      V = V(:,idx);
      D = D(idx,idx);
    end
    Dvals = diag(D).';
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

  nuclearHeader = 'ELECTRIC AND MAGNETIC HYPERFINE STRUCTURE';
  HyperfineEFGData = false;
  for k = kstart:kend
    if strcmp(L{k},nuclearHeader)
      HyperfineEFGData = true;
      break
    end
  end

  if HyperfineEFGData
    iAtom = 0;
    hfc = cell(1,nAtoms);
    efg = cell(1,nAtoms);
    while k<=kend
      if regexp(L{k},'^\s*Nucleus\s*')
        iAtom = sscanf(L{k}(9:end),'%d',1)+1;
      elseif regexp(L{k},'^\s*Raw HFC matrix\s*')
        idx = k+1;
        if strncmp(L{k+1}(2:4),'---',3)
          idx = k+2;
        end
        hfc_(1,:) = sscanf(L{idx},'%f %f %f');
        hfc_(2,:) = sscanf(L{idx+1},'%f %f %f').';
        hfc_(3,:) = sscanf(L{idx+2},'%f %f %f').';
        hfc{iAtom} = hfc_;
        k = k+3;
      elseif regexp(L{k},'^\s*Raw EFG matrix\s*')
        idx = k+1;
        efg_(1,:) = sscanf(L{idx+1},'%f %f %f');
        efg_(2,:) = sscanf(L{idx+2},'%f %f %f').';
        efg_(3,:) = sscanf(L{idx+3},'%f %f %f').';
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
        el = Element{iAtom};
        Qmom = nucqmom('17O')*barn;
        I = nucspin('17O');
        e2qQh = echarge*Qmom*eq(3)/planck/1e6;
        K = e2qQh/(4*I)/(2*I-1);
        Q{iAtom} = K*[1-eta,1+eta,-2];
        QFrame{iAtom} = eulang(R);
        k = k+3;
      end
      k = k+1;
    end
  end

  data(iStructure).hfc = hfc;
  data(iStructure).efg = efg;

end % for iStructure = 1:nStructures

% Copy relevant data to spin system structure
Sys = data;

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
