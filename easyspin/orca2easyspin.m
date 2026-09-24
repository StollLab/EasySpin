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
%              In addition to the spin Hamiltonian parameters (S, g, D, A,
%              Q, frames, Nucs, NucsIdx), it contains the total charge
%              (Sys.charge), the element symbols of all atoms
%              (Sys.Elements), and the atom coordinates in Angstrom (Sys.xyz).
%              These are included if they are present in the file.
%    Sys.data  contains additional data read from the output file
%              (coordinates, charge, electric field gradients, etc)
%
%  Besides the main text-formatted output file, ORCA also generates an
%  additional file that contains atomic coordinates and calculated
%  properties such as g and A matrices, Q tensors, etc. This property file
%  is text-based and ends in .property.txt (ORCA 6 and later). Before
%  ORCA 5, the property file was binary and had extension .prop.
%  orca2easyspin can read the main output file or either of these property
%  files. The text-based property files from ORCA 5 (_property.txt) are
%  not supported.
%
%  Examples:
%    Sys = orca2easyspin('nitroxide.out')   % all ORCA versions
%    Sys = orca2easyspin('nitroxide.property.txt')   % ORCA v6 and later
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

if output_ext==".prop"
  % binary property file (ORCA versions < 5)
  readmode = 'propbin';
  fileType = 'binary ORCA property';
elseif output_ext==".txt" && endsWith(output_name,".property")
  % text-based property file (ORCA versions >= 6)
  readmode = 'proptxtv6';
  fileType = 'text-based ORCA property';
elseif output_ext==".txt" && endsWith(output_name,"_property")
  error('ORCA 5 property files (_property.txt) are not supported. Provide the main output file instead.');
else
  % main ORCA output file
  readmode = 'mainout';
  fileType = 'ORCA output';
end

if ~exist(OrcaOutput,'file')
  error('Cannot access %s file %s.',fileType,OrcaOutput);
end


% Read properties from main or property output files
%--------------------------------------------------------------------------
switch readmode
  case "mainout"
    Sys = orca2easyspin_mainout(OrcaOutput);
  case "propbin"
    Sys = orca2easyspin_propbin(OrcaOutput);
  case "proptxtv6"
    Sys = orca2easyspin_proptxtv6(OrcaOutput);
end

% Apply hyperfine cutoff
%--------------------------------------------------------------------------
Sys = nucspinhftrim(Sys,HyperfineCutoff);

end
%==========================================================================


% Remove all nuclei with hyperfine coupling strength below a threshold.
% Nuclei without hyperfine data (all-zero A) but with quadrupole data are kept.
function Sys = nucspinhftrim(Sys,HyperfineCutoff)
if ~isfield(Sys,'Nucs')
  return
end
perNucFields = {'A','AFrame','Q','QFrame'};
for iSys = 1:numel(Sys)
  if isempty(Sys(iSys).Nucs)
    continue
  end
  Nucs = nucstring2list(Sys(iSys).Nucs);
  nNucs = numel(Nucs);
  if isfield(Sys,'A') && ~isempty(Sys(iSys).A)
    Amax = max(abs(Sys(iSys).A),[],2).';
  else
    Amax = zeros(1,nNucs);
  end
  if isfield(Sys,'Q') && ~isempty(Sys(iSys).Q)
    hasQ = any(Sys(iSys).Q,2).';
  else
    hasQ = false(1,nNucs);
  end
  keep = Amax>abs(HyperfineCutoff) | (Amax==0 & hasQ);
  Sys(iSys).Nucs = nuclist2string(Nucs(keep));
  Sys(iSys).NucsIdx = Sys(iSys).NucsIdx(keep);
  for f = perNucFields
    if isfield(Sys,f{1}) && ~isempty(Sys(iSys).(f{1}))
      Sys(iSys).(f{1}) = Sys(iSys).(f{1})(keep,:);
    end
  end
end
end
