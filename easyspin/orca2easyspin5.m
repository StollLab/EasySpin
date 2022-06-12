function Sys = orca2easyspin5(OrcaFileName)

if nargin<1
  %OrcaFileName = 'hydroxyl_HO_property.txt';
  OrcaFileName = 'hydroxyl_Honly_property.txt';
end

% Read entire file
%-------------------------------------------------------------------------------
fh = fopen(OrcaFileName);
allLines = textscan(fh,'%s','whitespace','','delimiter',newline);
allLines = allLines{1};
fclose(fh);

% Find all headers
%-------------------------------------------------------------------------------
headers = regexp(allLines,'-\s!(.*)!\s-','tokens');
headeridx = find(~cellfun(@isempty,headers));
headers = cellfun(@(x) x{1}{1},headers(headeridx),'UniformOutput',false);

% Find PROPERTIES sections
%-------------------------------------------------------------------------------
sections = regexp(allLines,'^\$\s(.*)','tokens');
sectionidx = find(~cellfun(@isempty,sections));
sections = cellfun(@(x) x{1}{1},sections(sectionidx),'UniformOutput',false);

% Parse PROPERTIES sections
%-------------------------------------------------------------------------------
for s = 1:numel(sectionidx)
  sidx = sectionidx(s);
  h = find(headeridx<sidx,1,'last');
  %description = sscanf(allLines{sidx+1}(17:end),'%d');
  geomidx = sscanf(allLines{sidx+2}(17:end),'%d');
  propidx = sscanf(allLines{sidx+3}(17:end),'%d');
  p = struct;
  switch sections{s}
    case 'Calculation_Info'
      startidx = 45;
      parsenum = @(sh,form)sscanf(allLines{sidx+sh}(startidx:end),form);
      p.multiplicity = parsenum(4,'%d');
      p.charge = parsenum(5,'%d');
      p.natoms = parsenum(6,'%d');
      p.nelectrons = parsenum(7,'%d');
      p.nfrozenelectrons = parsenum(8,'%d');
      p.ncorrelectrons = parsenum(9,'%d');
      p.nbasis = parsenum(10,'%d');

    case 'SCF_Energy'
      p.E = sscanf(allLines{sidx+4}(21:end),'%f');

    case 'EPRNMR_GTensor'
      % Bug in ORCA 5.0.1: raw g tensor is not in property file.
      % So we reconstruct (symmetrized) g tensor matrix from eigenvectors and
      % eigenvales.
      %g = read3by3array(allLines(idx0+(3:5)));
      R = read3by3array(allLines(sidx+(8:10)));
      gpv = sscanf(allLines{sidx+13},'%f').';
      gpv = gpv(2:4);
      g = R*diag(gpv)*R.';
      p.g = g;

    case 'EPRNMR_ATensor'
      st = find(allLines{sidx+4}==' ',1,'last');
      nItems = sscanf(allLines{sidx+4}(st:end),'%d');
      for iNuc = 1:nItems
        idx0 = sidx + 7 + (iNuc-1)*18;
        A = read3by3array(allLines(idx0+(6:8)));
        [atomID,element] = strtok(allLines{idx0}(10:end));
        atomID = sscanf(atomID,'%d');
        element = strtrim(element);
        isotope = sscanf(allLines{idx0+1}(10:end),'%d');
        p.Nucs{atomID+1} = sprintf('%d%s',isotope,element);
        p.A{atomID+1} = A;
        [R,Apv] = eig(A);
        p.Apv{atomID+1} = diag(Apv).';
        p.AFrame{atomID+1} = R;
      end

    case 'EPRNMR_QTensor'
      st = find(allLines{sidx+4}==' ',1,'last');
      nItems = sscanf(allLines{sidx+4}(st:end),'%d');
      for iNuc = 1:nItems
        idx0 = sidx + 8 + (iNuc-1)*19;
        % Bugs in ORCA v5.0.1:
        %  - Stored hyperfine tensor instead of Q tensor
        %  - Eigenvalues are of EFG tensor, not Q tensor
        %Q = read3by3array(allLines(idx0+(6:8)));
        R = read3by3array(allLines(idx0+(11:13)));
        efgpv = sscanf(allLines{idx0+16},'%f').';
        efgpv = efgpv(2:4); % in atomic units

        [atomID,element] = strtok(allLines{idx0}(10:end));
        atomID = sscanf(atomID,'%d');
        element = strtrim(element);
        isotope = sscanf(allLines{idx0+1}(10:end),'%d');
        p.Nucs{atomID+1} = sprintf('%d%s',isotope,element);
        Qpv = efg2Q(efgpv,element,'au');
        Q = R*diag(Qpv)*R.';
        p.Qpv{atomID+1} = Qpv;
        p.QFrame{atomID+1} = R;
        p.Q{atomID+1} = Q;
      end

    case 'EPRNMR_DTensor'
      error('Import of D tensor from EPRNMR_DTensor section is not implemented yet.');

  end

  if strcmp(headers{h},'PROPERTIES')
    PROPERTIES{geomidx,propidx}.(sections{s}) = p;
  end

end

if numel(PROPERTIES)==1
  PROPERTIES = PROPERTIES{1,1};
else
  error('Cannot currenly handle siutations with multiple geometries.');
end

PROPERTIES.EPRNMR_ATensor.A
PROPERTIES.Calculation_Info.natoms

% Build EasySpin spin system
%-------------------------------------------------------------------------------
Sys.S = (PROPERTIES.Calculation_Info.multiplicity-1)/2;
Sys.g = PROPERTIES.EPRNMR_GTensor.g;

nAtoms = PROPERTIES.Calculation_Info.natoms;
if isfield(PROPERTIES,'EPRNMR_ATensor')
  Apv = PROPERTIES.EPRNMR_ATensor.Apv;
  AFrame = PROPERTIES.EPRNMR_ATensor.AFrame;
  if numel(Apv)<nAtoms, Apv{nAtoms} = []; end
  if numel(AFrame)<nAtoms, AFrame{nAtoms} = []; end
else
  Apv = cell(1,nAtoms);
  AFrame = cell(1,nAtoms);
end
if isfield(PROPERTIES,'EPRNMR_QTensor')
  Qpv = PROPERTIES.EPRNMR_QTensor.Qpv;
  QFrame = PROPERTIES.EPRNMR_QTensor.QFrame;
  if numel(Qpv)<nAtoms, Qpv{nAtoms} = []; end
  if numel(QFrame)<nAtoms, QFrame{nAtoms} = []; end
else
  Qpv = cell(1,nAtoms);
  QFrame = cell(1,nAtoms);
end

idx = 0;
for iAtom = 1:PROPERTIES.Calculation_Info.natoms
  if ~isempty(Apv{iAtom}) || ~isempty(Qpv{iAtom})
    idx = idx + 1;
  else
    continue
  end
  if ~isempty(Apv{iAtom}) 
    Sys.A(idx,:) = Apv{iAtom};
    Sys.AFrame(idx,:) = eulang(AFrame{iAtom});
  else
    Sys.A(idx,:) = [NaN NaN NaN];
  end

  if ~isempty(Qpv{iAtom})
    Sys.Q(idx,:) = Qpv{iAtom};
    Sys.QFrame(idx,:) = eulang(QFrame{iAtom});
  else
    Sys.Q(idx,:) = [NaN NaN NaN];
  end
end

end
%===============================================================================


function T = read3by3array(lines)
T1 = sscanf(lines{1},'%f');
T2 = sscanf(lines{2},'%f');
T3 = sscanf(lines{3},'%f');
T = [T1 T2 T3];
T = T(2:end,:).';
end


function Qpv = efg2Q(efgpv,elem,units)
% Given the principal values of the electric field gradient (EFG) tensor,
% calculate the principal values for the nuclear quadrupole (Q) tensor for
% the most naturally abundant isotope with spin > 1/2.
%
% Input:
%   efgpv   principal values of EFG tensor, in SI units (V/m^2) or atomic
%           units (Eh/e/a0^2)
%   elem    element number
%   units   unit of EFG tensor principal values, 'SI' or 'au'
%
% Output:
%   Qpv     principal values of Q tensor, in MHz

switch units
  case 'SI'
    % nop
  case 'au'
    efgpv = efgpv*(hartree/evolt/bohrrad^2);  % Eh/e/a0^2 -> V/m^2
  otherwise
    error('Unrecognized unit ''%s'' for EFG. Use ''SI'' or ''au''.',unit);
end

% List of quadrupole reference isotopes for all elements
%  (most naturally abundant with I>1/2)
qrefMassNo = [...
  2  0  7  9  11 0  14 17 0  21 ... % 1-10
  23 25 27 0  0  33 35 0  39 43 ... % 11-20
  45 47 51 53 55 0  59 61 63 67 ...% 21-30
  69 73 75 0  79 83 85 87 0  91 ... % 31-40
  93 95 0  101 0 105 0 0 115 0 ... % 41-50
  121 0 127 131 133 137 139 0 141 143 ... % 51-60
  0 147 153 157 159 163 165 167 0 173 ... % 61-70
  175 177 181 0 187 189 193 0 197 201 ... % 71-80
  0 0 209 0 0 0 0 0 227 0 ... % 81-90
  0 235 237 0 243 0 0 0 0 0 ... % 91-100
  0 0 0 0 0 0 0 0 0 0 ... % 101-110
  0 0 0 0 0 0 0 0]; % 111-118
massNo = qrefMassNo(elem);

if massNo==0
  Qpv = [0 0 0];
else
  qrefIso = sprintf('%d%s',massNo,elementno2symbol(elem));
  eQ = echarge*nucqmom(qrefIso)*barn; % nuclear electric quadrupole moment, SI unit (C m^2)
  I = nucspin(qrefIso);
  Qpv_ = eQ/2/I/(2*I-1)*efgpv; % quadrupole tensor principal values, SI unit (J)
  Qpv = Qpv_/planck/1e6; % J -> MHz
end

end
