% orca2easyspin   Import magnetic spin Hamiltonian parameters from ORCA
%
%  Sys = orca2easyspin(OrcaFileName)
%  Sys = orca2easyspin(OrcaFileName,HyperfineCutoff)
%
%  Loads the spin Hamiltonian parameters from the binary ORCA property
%  output file given in OrcaFileName and returns them as an EasySpin
%  spin system structure Sys.
%
%  Input:
%    OrcaFileName     file name of the main ORCA output file
%    HyperfineCutoff  cutoff for hyperfine coupling (MHz)
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
%    Sys = orca2easyspin('nitroxide.out')   % all versions
%    Sys = orca2easyspin('nitroxide.prop')   % before ORCA 5
%    Sys = orca2easyspin('nitroxide_property.txt')   % ORCA 5 and later
%    Sys = orca2easyspin('nitroxide.out',0.5)  % 0.5 MHz hf cutoff

function Sys = orca2easyspin(OrcaOutput,HyperfineCutoff)

if nargin==0 && nargout==0
  help(mfilename);
  return
end

if nargin<2
  HyperfineCutoff = 0;  % MHz
end

% Detect type of ORCA file provided
%--------------------------------------------------------------------------
[output_path,output_name,output_ext] = fileparts(OrcaOutput);

if output_ext==".prop"
  % binary properry file (ORCA versions < 5.0)
  readmode = 'propbin';
  mainOutputFile = fullfile(output_path,output_name);
  binaryPropFile = fullfile(output_path,[output_name output_ext]);
  textPropFile = '';
elseif output_ext==".txt" && numel(output_name)>9 && output_name(end-8:end)=="_property"
  % text property file (ORCA versions >= 5.0)
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
existBin = exist(binaryPropFile,'file');
existTxt = exist(textPropFile,'file');

if readmode=="mainout" && ~existMain
  error('Cannot access ORCA output file %s.',mainOutputFile);
end
if readmode=="propbin" && ~existBin
  error('Cannot access ORCA property file %s.',binaryPropFile);
end
if readmode=="`proptxt" && ~existTxt
  error('Cannot access ORCA property file %s.',textPropFile);
end

% Determine ORCA version, if possible
%--------------------------------------------------------------------------
if existMain
  vOrca = getOrcaVersion(mainOutputFile);

  if (vOrca=="5.0.1" || vOrca=="5.0.0") && readmode=="proptxt"
    error('Cannot read text-based property file for ORCA versions 5.0.0 and 5.0.1, since it''s buggy. Use main file instead.');
  end
end


% Read properties from main or property output files
%--------------------------------------------------------------------------
switch readmode
  case "mainout"
    Sys = orca2easyspin_maintxt(mainOutputFile,HyperfineCutoff);
  case "propbin"
    Sys = orca2easyspin_propbin(binaryPropFile,HyperfineCutoff);
  case "proptxt"
    Sys = orca2easyspin_proptxt(textPropFile,HyperfineCutoff);
end

end
%==========================================================================

% Determine ORCA version by looking through the top of the text-based
% output file for a line containing "Program Version x.y.z"
function OrcaVersion = getOrcaVersion(mainOutputFile)
OrcaVersion = '';
maxLines = 50; % limit search to initial lines
if ~isempty(mainOutputFile)
  f = fopen(mainOutputFile);
  idx = 0;
  while isempty(OrcaVersion) && idx<maxLines && ~feof(f)
    idx = idx+1;
    thisLine = fgetl(f);
    OrcaVersion = regexp(thisLine,'\d+\.\d+\.\d+','match','once');
  end
  fclose(f);
end
end
