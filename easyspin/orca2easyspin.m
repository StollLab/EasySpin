% orca2easyspin   Import magnetic spin Hamiltonian parameters from ORCA
%
%  Sys = orca2easyspin(OrcaFileName)
%  Sys = orca2easyspin(OrcaFileName,HyperfineCutoff)
%
%  Loads the spin Hamiltonian parameters from the binary ORCA property
%  output file given in OrcaFileName and returns them as an EasySpin
%  spin system structure Sys.
%
%  Besides the main text-formatted output file, ORCA also generates a
%  binary .prop file that contains atomic coordinates and calculated
%  properties such as g and A matrices, Q tensors, etc. orca2easyspin
%  reads in this file, and not the text-formatted output file.
%
%  If HyperfineCutoff (a single value, in MHz) is given, all
%  nuclei with hyperfine coupling equal or smaller than that
%  value are omitted from the spin system. If not given, it
%  is set to zero, and all nuclei with non-zero hyperfine
%  coupling are included.
%
%  Examples:
%    Sys = orca2easyspin('nitroxide.prop')
%    Sys = orca2easyspin('nitroxide.prop',0.5)  % 0.5 MHz hf cutoff

function Sys = orca2easyspin(propfilename,HyperfineCutoff)

if nargin==0 && nargout==0, help(mfilename); return; end

if nargin==1, HyperfineCutoff = 0; end

% Whether or not to use least-squares fitting for
% rotation matrix -> Euler angle conversion using eulang()
skipFitting = true;

DebugMode = false;

Sys = [];

PropertyIDtype = 'int32'; % datatype for property ID in .prop file
rowcoltype = 'int32';
datatype = 'float64';
terminateID = -1; % property ID value that indicates end of file

% Property IDs have to be between 0 and 53 (as of ORCA 3.0.2)
maxPropertyID = 53;

% Get base name of filename and compose .prop filename
%-------------------------------------------------------------------------------
[pf_path,pf_name,~] = fileparts(propfilename);
propfilename = fullfile(pf_path,[pf_name '.prop']);

% Determine correct machine format (little-endian vs. big-endian)
%-------------------------------------------------------------------------------
MachineFormat = 'l'; % little-endian
[f,errmsg] = fopen(propfilename,'r',MachineFormat);
if f<0
  error('Could not open file %s: %s',propfilename,errmsg);
end
firstPropertyID = fread(f,1,PropertyIDtype);
fclose(f);

% quit if property file is empty (4 bytes, all 0xFF)
if firstPropertyID==terminateID
  return
end

if firstPropertyID>maxPropertyID
  % Clearly, if the read ID is outside the valid range, the original
  % MachineFormat guess was wrong. Try the other one.
  MachineFormat = 'b'; % big-endian
end


% Read in all EPR-relevant properties
%-------------------------------------------------------------------------------
f = fopen(propfilename,'r',MachineFormat);

S = [];
xyz = [];
Charge = [];
Dpv = [];
DFrame = [];
gpv = [];
gFrame = [];
Apv = [];
AFrame = [];
efg = [];
efgFrame = [];
rho0 = [];
Atoms = [];

while ~feof(f)
  
  PropertyID = fread(f,1,PropertyIDtype);
  
  if DebugMode
    fprintf('%2d - %s\n',PropertyID,propertystring(PropertyID)); %#ok<UNRCH>
  end
  
  if PropertyID==terminateID, break; end
  
  if PropertyID<0 || PropertyID>maxPropertyID
    error('Unknown property encountered (PropertyID = %d)',PropertyID);
  end
  
  nRows = fread(f,1,rowcoltype);
  nColumns = fread(f,1,rowcoltype);
  nElements = nRows*nColumns;
  
  data = fread(f,nElements,datatype);
  data = reshape(data,nColumns,nRows);

  switch PropertyID

    % Charge -------------------------------------------
    case 40
      Charge = data(1);
      
    % Spin multiplicity --------------------------------
    case 41
      S = (data(1)-1)/2;
  
    % Coordinates --------------------------------------
    case {17,42} % 17 up to ORCA 3.0.2, then 42
      % only included if HF are calculated
      xyz = reshape(data,3,[]).'; % in Bohr radii
      xyz = xyz*bohrrad/1e-10;    % conversion to Angstrom
    
    % Element numbers ----------------------------------
    case {18,43} % 18 up to ORCA 3.0.2, then 43
      Atoms = data.';
      
    % g matrix -----------------------------------------
    case {5, 14, 23, 30, 27}
      gpv = data(1:3).';
      if all(gpv==0)
        warning('Property entry with ID %d (%s) is all zero.',PropertyID,propertystring(PropertyID));
        gFrame = [0 0 0];
      else
        R = reshape(data(4:12),3,3);
        %g = R*diag(gpv)*R.';
        gFrame = eulang(R.',skipFitting);
      end
      
    % D tensor -----------------------------------------
    case {4, 13, 22, 29, 36}
      %Dval = data(1);
      %EoverD = data(2);
      Dpv = data(3:5).';       % cm^-1
      Dpv = Dpv*1e-4*clight;   % cm^-1 -> MHz
      R = reshape(data(6:14),3,3);
      %D = R*diag(Dpv)*R.';
      DFrame = eulang(R.',skipFitting);
      
    % A matrices ---------------------------------------
    case {6, 15, 24, 31, 38}
      nucIdx = data(1,:)+1;
      %aiso = data(2,:);
      for iNuc=1:numel(nucIdx)
        idx = nucIdx(iNuc);
        Apv(idx,1:3) = data(3:5,iNuc).'; %#ok<AGROW>
        R = reshape(data(6:14,iNuc),3,3);
        %A = R*diag(Apv)*R.';
        AFrame(idx,1:3) = eulang(R.',skipFitting).'; %#ok<AGROW>
      end
    
    % EFG tensors ----------------------------------------
    case {7, 16, 25, 32, 39}
      nucIdx = data(1,:)+1;
      for iNuc = numel(nucIdx):-1:1
        idx = nucIdx(iNuc);
        efg_au = data(2:4,iNuc).'; % EFG, atomic unit (Eh/e/a0^2)
        efg(idx,1:3) = efg_au*(hartree/echarge/bohrrad^2); % -> SI unit (V/m^2)
        R = reshape(data(5:13,iNuc),3,3).';
        efgFrame(idx,1:3) = eulang(R.',skipFitting);
      end
    
    % Spin densities at nuclei -------------------------
    case 9
      nucIdx = data(1,:) + 1; % ORCA is 0-based, MATLAB is 1-based
      rho0_ = data(2,:); % atomic units (a0^-3)
      rho0(nucIdx) = rho0_*(bohrrad/1e-10)^3; %#ok<AGROW> % conversion to angstrom^-3 
      
    otherwise
      % skip other properties (dipole moment, polarizability, etc)
  end
end
fclose(f);

nAtoms = numel(Atoms);

% Compile spin system
%-------------------------------------------------------------------------------
if isempty(S)
  % spin is not provided by the prop file
else
  Sys.S = S;
end
if ~isempty(xyz)
  Sys.xyz = xyz;
  Sys.Atoms = Atoms;
end
if ~isempty(Charge)
  Sys.Charge = Charge;
end
if ~isempty(gpv)
  Sys.g = gpv;
  Sys.gFrame = gFrame;
end
if ~isempty(Dpv)
  Sys.D = Dpv;
  Sys.DFrame = DFrame;
end

% Pad with zeros if necessary
anyHyperfine = ~isempty(Apv);
anyQuadrupole = ~isempty(efg);
if anyHyperfine
  if size(Apv,1)<nAtoms, Apv(nAtoms,:) = 0; end
  if size(AFrame,1)<nAtoms, AFrame(nAtoms,:) = 0; end
end
if anyQuadrupole
  if size(efg,1)<nAtoms, efg(nAtoms,:) = 0; end
  if size(efgFrame,1)<nAtoms, efgFrame(nAtoms,:) = 0; end
  Qpv = zeros(size(efg));
end

if nAtoms>0
  
  % Convert electric field gradient to Q tensor principal values
  if anyQuadrupole
    for iAtom = nAtoms:-1:1
      if ~any(efg(iAtom,:)), continue; end
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
      massNo = qrefMassNo(Atoms(iAtom));
      if massNo==0
        Qpv(iAtom,:) = [0 0 0];
      else
        qrefIso = sprintf('%d%s',massNo,elementno2symbol(Atoms(iAtom)));
        eQ = echarge*nucqmom(qrefIso)*barn; % nuclear electric quadrupole moment, SI unit (C m^2)
        I = nucspin(qrefIso);
        Qpv_ = eQ/2/I/(2*I-1)*efg(iAtom,:); % quadrupole tensor principal values, SI unit (J)
        Qpv(iAtom,:) = Qpv_/planck/1e6; % J -> MHz
      end
    end
  end
  
  % Determine which nuclei are above hyperfine cutoff
  hfkeep = false(1,nAtoms);
  if anyHyperfine
    Amax = max(abs(Apv),[],2);
    hfkeep = Amax>HyperfineCutoff;
  end
  
  % Build Sys.Nucs
  NucStr = [];
  keepidx = 1:nAtoms;
  for iAtom = 1:nAtoms
    if ~hfkeep(iAtom)
        keepidx(iAtom) = false;
        continue
    end
    NucStr = [NucStr ',' elementno2symbol(Atoms(iAtom))]; %#ok<AGROW>
  end
  if ~isempty(NucStr)
    NucStr(1) = [];
  end
  Sys.Nucs = NucStr;
  Sys.NucsIdx = find(keepidx);
  
  % Build Sys.A and Sys.AFrame
  if anyHyperfine
    Sys.A = Apv(hfkeep,:);
    Sys.AFrame = AFrame(hfkeep,:);
  end
  
  % Build Sys.Q and Sys.QFrame
  if anyQuadrupole
    Sys.Q = Qpv(hfkeep,:);
    Sys.QFrame = efgFrame(hfkeep,:);
  end
  
end
%===============================================================================


% Based on char *PropertyID definition in qcpropfl.h from ORCA source code
% (unclear which version...)
function s = propertystring(ID)
switch ID
  case  0, s = 'SCF_DIP';
  case  1, s = 'SCF_POL';
  case  2, s = 'SCF_ABS';
  case  3, s = 'SCF_CD';
  case  4, s = 'SCF_D';
  case  5, s = 'SCF_G';
  case  6, s = 'SCF_A';
  case  7, s = 'SCF_Q';
  case  8, s = 'SCF_SHIFT';
  case  9, s = 'SCF_RHO_0';
  
  case 10, s = 'MRCI_DIP';
  case 11, s = 'MRCI_ABS';
  case 12, s = 'MRCI_CD';
  case 13, s = 'MRCI_D';
  case 14, s = 'MRCI_G';
  case 15, s = 'MRCI_A';
  case 16, s = 'MRCI_Q';
  
  case 17, s = 'XYZ';
  case 18, s = 'ATNO';
  
  case 19, s = 'CASSCF_DIP';
  case 20, s = 'CASSCF_ABS';
  case 21, s = 'CASSCF_CD';
  case 22, s = 'CASSCF_D';
  case 23, s = 'CASSCF_G';
  case 24, s = 'CASSCF_A';
  case 25, s = 'CASSCF_Q';
  
  case 26, s = 'NEVPT2_DIP';
  case 27, s = 'NEVPT2_ABS';
  case 28, s = 'NEVPT2_CD';
  case 29, s = 'NEVPT2_D';
  case 30, s = 'NEVPT2_G';
  case 31, s = 'NEVPT2_A';
  case 32, s = 'NEVPT2_Q';
  
  case 33, s = 'CASSCFCUSTOME_DIP';
  case 34, s = 'CASSCFCUSTOME_ABS';
  case 35, s = 'CASSCFCUSTOME_CD';
  case 36, s = 'CASSCFCUSTOME_D';
  case 37, s = 'CASSCFCUSTOME_G';
  case 38, s = 'CASSCFCUSTOME_A';
  case 39, s = 'CASSCFCUSTOME_Q';

  case 40, s = 'QDNEVPT2_DIP';
  case 41, s = 'QDNEVPT2_ABS';
  case 42, s = 'QDNEVPT2_CD';
  case 43, s = 'QDNEVPT2_D';
  case 44, s = 'QDNEVPT2_G';
  case 45, s = 'QDNEVPT2_A';
  case 46, s = 'QDNEVPT2_Q';
  
  case -1, s = 'termination';
  
  otherwise, s = 'unknown';
  
end
