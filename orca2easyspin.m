% orca2easyspin   Import magnetic parameters from ORCA calculation
%
%  Sys = orca2easyspin(OrcaFileName)
%  Sys = orca2easyspin(OrcaFileName,HyperfineCutoff)
%
%  Loads the magnetic properties from the binary ORCA property output
%  file given in OrcaFileName and returns them as an EasySpin spin
%  system structure Sys.
%
%  Besides the text-formatted output file, ORCA also generates a
%  binary .prop file that contains atomic coordinates and calculated
%  properties such as g  and A matrices, Q tensors, etc. orca2easyspin
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

if (nargin==0)&&(nargout==0), help(mfilename); return; end

if (nargin==1), HyperfineCutoff = 0; end

Sys = [];

PropertyIDtype = 'int32'; % datatype for property ID in .prop file
rowcoltype = 'int32';
datatype = 'float64';
terminateID = -1; % property ID value that indicates end of file

% Get base name of filename and compose .prop filename
%----------------------------------------------------------------
[pf_path,pf_name,pf_ext] = fileparts(propfilename);
propfilename = fullfile(pf_path,[pf_name '.prop']);

% Determine correct machine format (little-endian vs. big-endian)
%----------------------------------------------------------------
MachineFormat = 'l'; % little-endian
[f,errmsg] = fopen(propfilename,'r',MachineFormat);
if (f<0)
  error('Could not open file %s: %s',propfilename,errmsg);
end
firstPropertyID = fread(f,1,PropertyIDtype);
fclose(f);

% quit if property file is empty (4 bytes, all 0xFF)
if (firstPropertyID==terminateID)
  return
end

if (firstPropertyID>100)
  % Property IDs have to be between 0 and 53 (as of ORCA 3.0.2)
  % Clearly, if the first one is larger, the original MachineFormat
  % guess was wrong. Try the other one.
  MachineFormat = 'b'; % big-endian
end

% Whether or not to use least-squares fitting for
% rotation matrix -> Euler angle conversion using eulang()
skipFitting = true;

% Read in all EPR-relevant properties
%--------------------------------------------------------------
f = fopen(propfilename,'r',MachineFormat);

S = [];
xyz = [];
Charge = [];
Dpv = [];
gpv = [];
Apv = [];
Qpv = [];
rho0 = [];

while ~feof(f)
  
  PropertyID = fread(f,1,PropertyIDtype);
  if (PropertyID==terminateID), break; end
  
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
      S = data(1);
  
    % Coordinates --------------------------------------
    case {17,42} % 17 up to ORCA 3.0.2, then 42
      xyz = reshape(data,3,[]).'; % in Bohr radii
      xyz = xyz*bohrrad/1e-10;    % conversion to Angstrom
    
    % Element numbers ----------------------------------
    case {18,43} % 18 up to ORCA 3.0.2, then 43
      Atoms = data.';
      
    % g matrix -----------------------------------------
    case {5, 14, 23, 30, 27}
      gpv = data(1:3).';
      R = reshape(data(4:12),3,3);
      %g = R*diag(gpv)*R.';
      gFrame = eulang(R.',skipFitting);
      
    % D tensor -----------------------------------------
    case {4, 13, 22, 29, 36}
      %Dval = data(1);
      %EoverD = data(2);
      Dpv = data(3:5);
      R = reshape(data(6:14),3,3);
      %D = R*diag(Dpv)*R.';
      DFrame = eulang(R.',skipFitting);
    
    % A matrices ---------------------------------------
    case {6, 15, 24, 31, 38}
      AnucIdx = data(1,:)+1;
      %aiso = data(2,:);
      for iNuc=1:numel(AnucIdx)
        idx = AnucIdx(iNuc);
        Apv(idx,1:3) = data(3:5,iNuc).';
        R = reshape(data(6:14),3,3);
        %A = R*diag(Apv)*R.';
        AFrame(idx,1:3) = eulang(R.',skipFitting).';
      end
    
    % Q tensors ----------------------------------------
    case {7, 16, 25, 32, 39}
      QnucIdx = data(1,:)+1;
      for iNuc = 1:numel(QnucIdx)
        idx = QnucIdx(iNuc);
        Qpv(idx,1:3) = data(2:4,iNuc).';
        R = reshape(data(5:13,iNuc),3,3).';
        %Q = R*diag(Qpv)*R.';
        QFrame(idx,1:3) = eulang(R.',skipFitting);
      end
    
    % Spin densities at nuclei -------------------------
    case 9
      nucIdx = data(1,:) + 1; % ORCA is 0-based, MATLAB is 1-based
      rho0(nucIdx) = data(2,:); % atomic units (a0^-3)
      rho0(nucIdx) = rho0(nucIdx)*(bohrrad/1e-10)^3; % conversion to Angstrom^-3
      
    otherwise
      % skip other properties (dipole moment, polarizability, etc)
  end
end

nAtoms = numel(Atoms);
if size(Apv,1)<nAtoms, Apv(nAtoms,:) = 0; end
if size(AFrame,1)<nAtoms, AFrame(nAtoms,:) = 0; end
if size(Qpv,1)<nAtoms, Qpv(nAtoms,:) = 0; end
if size(QFrame,1)<nAtoms, QFrame(nAtoms,:) = 0; end

fclose(f);

% Compile spin system
%---------------------------------------------------------------
if isempty(S)
  % spin is not provided by the prop file
  Sys.S = 1/2;
  %fprintf('Spin not provided in the ORCA %s.prop file. Assuming S = 1/2.\n',pf_name);
else
  Sys.S = S;
end
Sys.xyz = xyz;
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

% Hyperfine cutoff
Amax = max(abs(Apv),[],2);
hfkeep = Amax>HyperfineCutoff;

% List of isotopes
NucStr = [];
for iAtom = 1:nAtoms
  if ~hfkeep(iAtom), continue; end
  NucStr = [NucStr ',' elementno2symbol(Atoms(iAtom))];
end
if ~isempty(NucStr)
  NucStr(1) = [];
end
Sys.Nucs = NucStr;

if ~isempty(Apv)
  Sys.A  = Apv(hfkeep,:);
  Sys.AFrame = AFrame(hfkeep,:);
end
if ~isempty(Qpv)
  Sys.Q  = Qpv(hfkeep,:);
  Sys.QFrame = QFrame(hfkeep,:);
end
