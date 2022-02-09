% orca2easyspin   Import magnetic spin Hamiltonian parameters from ORCA
%
%  Sys = orca2easyspin(OrcaFileName)
%  Sys = orca2easyspin(OrcaFileName,HyperfineCutoff)
%  Sys = orca2easyspin(OrcaFileName,HyperfineCutoff,readmode)
%
%  Loads the spin Hamiltonian parameters from the binary ORCA property
%  output file given in OrcaFileName and returns them as an EasySpin
%  spin system structure Sys.
%
%  Input:
%    OrcaFileName     file name of the main ORCA output file
%    HyperfineCutoff  cutoff for hyperfine coupling (MHz)
%    readmode         'main' - read main file (slow)
%                     'prop' - read property file (fast)
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
%    Sys = orca2easyspin('nitroxide.out')
%    Sys = orca2easyspin('nitroxide.out')
%    Sys = orca2easyspin('nitroxide.out',0.5)  % 0.5 MHz hf cutoff

function Sys = orca2easyspin(OrcaOutput,HyperfineCutoff,readmode)

if nargin==0 && nargout==0
  help(mfilename);
  return
end

if nargin<2
  HyperfineCutoff = 0;  % MHz
end

if nargin<3
  readmode = "prop";
end
readMainFile = readmode=="main";

% Deermine main ORCA output filename
%--------------------------------------------------------------------------
[output_path,output_name,output_ext] = fileparts(OrcaOutput);

if output_ext==".prop"
  mainOutputFile = '';
else
  mainOutputFile = fullfile(output_path,[output_name output_ext]);
end

% Determine ORCA version
%--------------------------------------------------------------------------
vOrca = getOrcaVersion(mainOutputFile);

if (vOrca=="5.0.1" || vOrca=="5.0.0") && ~readMainFile
  error('Cannot read text-based property file for ORCA versions 5.0.0 and 5.0.1, since it''s buggy.');
end


% Determine version-appropriate property file if needed
%--------------------------------------------------------------------------
textPropFile = fullfile(output_path,[output_name '_property.txt']);
if ~exist(textPropFile,'file')
  textPropFile = '';
end
binaryPropFile = fullfile(output_path,[output_name '.prop']);
if ~exist(binaryPropFile,'file')
  binaryPropFile = '';
end


% Read properties from main or property output files
%--------------------------------------------------------------------------
if readMainFile
  Sys = orca2easyspin_maintxt(mainOutputFile,HyperfineCutoff);
else
  if ~isempty(binaryPropFile)
    Sys = orca2easyspin_propbin(binaryPropFile,HyperfineCutoff);
  elseif ~isempty(textPropFile)
    Sys = orca2easyspin_proptxt(textPropFile,HyperfineCutoff);
  end
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
