% orca2easyspin   Import spin Hamiltonian parameters from ORCA
%
%  Sys = orca2easyspin(OrcaFileName)
%  Sys = orca2easyspin(OrcaFileName,HyperfineCutoff)
%
%  Loads the spin Hamiltonian parameters from the ORCA output file
%  given in OrcaFileName and returns them as an EasySpin spin system
%  structure Sys. If the output file contains multiple structures, Sys
%  is an array of spin system structures.
%
%  Input:
%    OrcaFileName     file name of the main ORCA output file
%    HyperfineCutoff  cutoff for hyperfine coupling (MHz)
%
%  Output:
%    Sys       spin system structure, or array of spin system structures
%    Sys.data  contains additional data read from the output file
%              (coordinates, charge, electric field gradients, etc)
%
%  Besides the main text-formatted output file, ORCA also generates an
%  additional file that contains atomic coordinates and calculated
%  properties such as g and A matrices, Q tensors, etc. This property file
%  is text-based and ends in _property.txt. Before ORCA 5, the property
%  file was binary and had extension .prop. orca2easyspin can read either
%  the main output file or the associated property file.
%
%  Examples:
%    Sys = orca2easyspin('nitroxide.out')   % all ORCA versions
%    Sys = orca2easyspin('nitroxide_property.txt')   % ORCA v5 and later
%    Sys = orca2easyspin('nitroxide.prop')   % ORCA prior to v5
%
%  If HyperfineCutoff (a single value, in MHz) is given, all nuclei with
%  hyperfine coupling equal or smaller than that value are omitted from
%  the spin system. If not given, it is set to zero, and all nuclei with
%  non-zero hyperfine coupling are included.
%
%  Example:
%    Sys = orca2easyspin('nitroxide.out',0.5)  % 0.5 MHz hyperfine cutoff

function Sys = orca2easyspin(OrcaOutput,HyperfineCutoff)

if nargin==0 && nargout==0
  help(mfilename);
  return
end

if nargin<2
  HyperfineCutoff = 0;  % MHz
end


% Detect type of ORCA output file provided
%--------------------------------------------------------------------------
[output_path,output_name,output_ext] = fileparts(OrcaOutput);

if output_ext==".txt" && numel(output_name)>9 && output_name(end-8:end)=="_property"
  error('orca2easyspin() currently does not support reading spin Hamiltonian parameters from _property.txt files. Provide the main output file instead.')
end

if output_ext==".prop"
  % binary property file (ORCA versions < 5)
  readmode = 'propbin';
  mainOutputFile = fullfile(output_path,output_name);
  binaryPropFile = fullfile(output_path,[output_name output_ext]);
  textPropFile = '';
elseif output_ext==".txt" && numel(output_name)>9 && output_name(end-8:end)=="_property"
  % text-based property file (ORCA versions >= 5)
  readmode = 'proptxt';
  mainOutputFile = fullfile(output_path,output_name(1:end-9));
  binaryPropFile = '';
  textPropFile = fullfile(output_path,[output_name output_ext]);
else
  % main ORCA output file
  readmode = 'mainout';
  mainOutputFile = fullfile(output_path,[output_name output_ext]);
  binaryPropFile = fullfile(output_path,[output_name output_ext '.prop']);
  textPropFile = fullfile(output_path,[output_name '_property.txt']);
end
existMain = exist(mainOutputFile,'file');
existPropBin = exist(binaryPropFile,'file');
existPropTxt = exist(textPropFile,'file');

if readmode=="mainout" && ~existMain
  error('Cannot access ORCA output file %s.',mainOutputFile);
end
if readmode=="propbin" && ~existPropBin
  error('Cannot access ORCA property file %s.',binaryPropFile);
end
if readmode=="proptxt" && ~existPropTxt
  error('Cannot access ORCA property file %s.',textPropFile);
end


% Block reading buggy property files (early ORCA 5 versions)
%--------------------------------------------------------------------------
if readmode=="proptxt"
  % Determine ORCA version, if possible
  vOrca = getOrcaVersion(mainOutputFile);
  buggyVersions = ["5.0.0", "5.0.1", "5.0.2", "5.0.3"];
  if any(vOrca==buggyVersions)
    error('Cannot read property file for ORCA version %s. Use main output file instead.',vOrca);
  end
end


% Read properties from main or property output files
%--------------------------------------------------------------------------
switch readmode
  case "mainout"
    Sys = orca2easyspin_maintxt(mainOutputFile);
  case "propbin"
    Sys = orca2easyspin_propbin(binaryPropFile);
  case "proptxt"
    Sys = orca2easyspin_proptxt(textPropFile);
end

% Apply hyperfine cutoff
%--------------------------------------------------------------------------
Sys = nucspinhftrim(Sys,HyperfineCutoff);

end
%==========================================================================


% Determine ORCA version by looking through the top of the text-based
% output file for a line containing "Program Version x.y.z"
function OrcaVersion = getOrcaVersion(mainOutputFile)
OrcaVersion = '';
if ~exist(mainOutputFile,'file')
  return
end
maxLines = 50; % limit search to initial lines
if ~isempty(mainOutputFile)
  fh = fopen(mainOutputFile);
  idx = 0;
  while isempty(OrcaVersion) && idx<maxLines && ~feof(fh)
    idx = idx + 1;
    thisLine = fgetl(fh);
    OrcaVersion = regexp(thisLine,'\d+\.\d+\.\d+','match','once');
  end
  fclose(fh);
end
end


% Remove all nuclei with hyperfine coupling strength below a threshold
function Sys = nucspinhftrim(Sys,HyperfineCutoff)
if isfield(Sys,'Nucs') && isfield(Sys,'A')
  for iSys = 1:numel(Sys)
    Amax = max(abs(Sys(iSys).A),[],2);
    keep = Amax > abs(HyperfineCutoff);
    if ~isfield(Sys,'Nucs') || isempty(Sys.Nucs)
      continue
    end
    Nucs = Sys(iSys).Nucs;
    if ischar(Nucs)
      Nucs = nucstring2list(Nucs);
    end
    Sys(iSys).NucsIdx = find(keep).';
    Sys(iSys).Nucs = nuclist2string(Nucs(keep));
    Sys(iSys).A = Sys(iSys).A(keep,:);
    if isfield(Sys,'AFrame')
      Sys(iSys).AFrame = Sys(iSys).AFrame(keep,:);
    end
    if isfield(Sys,'Q')
      Sys(iSys).Q = Sys(iSys).Q(keep,:);
    end
    if isfield(Sys,'QFrame')
      Sys(iSys).QFrame = Sys(iSys).QFrame(keep,:);
    end
  end
end
end
