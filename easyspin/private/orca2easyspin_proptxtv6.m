% orca2easyspin_proptxtv6   Read EPR properties from ORCA 6 property text file
%
%   [Sys,data] = orca2easyspin_proptxtv6(propfile)
%
%   Reads EPR-relevant information from a text-based ORCA property file
%   (*.property.txt, format introduced in ORCA 6) and constructs an EasySpin
%   spin system structure from it.
%
%   Input:
%     propfile   file name of the ORCA property file (*.property.txt)
%
%   Output:
%     Sys        spin system structure (array if the file contains several
%                geometries)
%                  .S        electron spin
%                  .charge   total charge
%                  .Elements element symbols of all atoms
%                  .xyz      atom coordinates (Å)
%                  .g, .gFrame, .D, .DFrame
%                  .Nucs, .NucsIdx, .A, .AFrame, .Q, .QFrame
%     data       structure with information as read from the file (one
%                element per geometry), in the units used in the file
%                  .OrcaVersion
%                  .Multiplicity, .Charge
%                  .Elements, .xyz (bohr)
%                  .g_matrix
%                  .d_raw (cm^-1)
%                  .ATensor.NucIdx, .Elem, .Isotope, .I,
%                  .PFAC (MHz), .ARaw (MHz)
%                  .EFGTensor.NucIdx, .Elem, .V (atomic unit, Eh/e/a0^2)
%                NucIdx are 1-based atom indices (0-based in the file).

function [Sys,data] = orca2easyspin_proptxtv6(propfile)

% File import and checks
%--------------------------------------------------------------------------
if ~exist(propfile,'file')
  error('Cannot find file %s.',propfile);
end
L = cellstr(readlines(propfile));
nLines = numel(L);

if any(contains(L(1:min(5,nLines)),'!PROPERTIES!'))
  error('%s is an ORCA property file in the pre-ORCA-6 format, which is not supported by this function.',propfile);
end

banner = '';
if nLines>=2
  banner = regexp(L{2},'^\*+\s*ORCA\s+(\d+\.\d+\.\d+)\s*\*+$','tokens','once');
end
if isempty(banner)
  error('%s is not an ORCA property file (*.property.txt).',propfile);
end

% Locate sections and assert proper $Title ... $End structure
isSectionLine = startsWith(L,'$');
sectionStart = find(isSectionLine & ~strcmp(L,'$End'));
sectionEnd = find(strcmp(L,'$End'));
if isempty(sectionStart) || numel(sectionStart)~=numel(sectionEnd) || ...
    any(sectionEnd<sectionStart) || any(sectionStart(2:end)<sectionEnd(1:end-1))
  error('%s is not a properly formed ORCA property file: section structure ($Title ... $End) is corrupt.',propfile);
end
sectionTitle = strtrim(extractAfter(L(sectionStart),1));

requiredSections = ["Calculation_Status","Geometry","Calculation_Info"];
for s = requiredSections
  if ~any(sectionTitle==s)
    error('%s is not a properly formed ORCA property file: $%s section is missing.',propfile,s);
  end
end

% Parse relevant sections
%--------------------------------------------------------------------------
relevantSections = ["Calculation_Status","Geometry","Calculation_Info",...
  "SCF_G_Tensor","SCF_D_Tensor","SCF_A_Tensor","SCF_EFG_Tensor"];

emptyData = struct('OrcaVersion',banner{1},'Multiplicity',[],'Charge',[],...
  'Elements',{{}},'xyz',[],'g_matrix',[],'d_raw',[],'ATensor',[],'EFGTensor',[]);
data = emptyData([]);

for iSection = 1:numel(sectionStart)
  title = sectionTitle{iSection};
  if ~any(title==relevantSections), continue; end

  P = parsesection(L(sectionStart(iSection)+1:sectionEnd(iSection)-1));
  if ~isfield(P,'GeometryIndex')
    error('Section $%s does not contain a geometry index.',title);
  end
  iGeom = P.GeometryIndex{1};
  if iGeom>numel(data)
    data(numel(data)+1:iGeom) = emptyData;
  end

  switch title
    case 'Calculation_Status'
      data(iGeom).OrcaVersion = getvalue(P,'version',title);

    case 'Geometry'
      [data(iGeom).Elements,data(iGeom).xyz] = getvalue(P,'CartesianCoordinates',title);

    case 'Calculation_Info'
      data(iGeom).Multiplicity = getvalue(P,'Mult',title);
      data(iGeom).Charge = getvalue(P,'Charge',title);

    case 'SCF_G_Tensor'
      data(iGeom).g_matrix = getvalue(P,'g_matrix',title);

    case 'SCF_D_Tensor'
      data(iGeom).d_raw = getvalue(P,'d_raw',title);

    case 'SCF_A_Tensor'
      A.NucIdx = getlist(P,'NUC',title) + 1;  % 0-based -> 1-based
      A.Elem = getlist(P,'Elem',title);
      A.Isotope = getlist(P,'Isotope',title);
      A.I = getlist(P,'I',title);
      A.PFAC = getlist(P,'PFAC',title);
      A.ARaw = cat(3,P.ARaw{:});
      checknuclearlists(A,{'Elem','Isotope','I','PFAC'},'ARaw',title);
      data(iGeom).ATensor = A;

    case 'SCF_EFG_Tensor'
      EFG.NucIdx = getlist(P,'NUC',title) + 1;  % 0-based -> 1-based
      EFG.Elem = getlist(P,'Elems',title);
      EFG.V = cat(3,P.V{:});
      checknuclearlists(EFG,{'Elem'},'V',title);
      data(iGeom).EFGTensor = EFG;
  end

end

for iGeom = 1:numel(data)
  if isempty(data(iGeom).xyz)
    error('No $Geometry section found for geometry %d.',iGeom);
  end
  if isempty(data(iGeom).Multiplicity)
    error('No $Calculation_Info section found for geometry %d.',iGeom);
  end
end

% Build spin system structures
%--------------------------------------------------------------------------
for iGeom = numel(data):-1:1
  Sys(iGeom) = buildspinsystem(data(iGeom));
end

end

%==========================================================================
function Sys = buildspinsystem(d)

Sys = struct;

% Spin, charge, geometry
Sys.S = (d.Multiplicity-1)/2;
Sys.charge = d.Charge;
Sys.Elements = d.Elements;
Sys.xyz = d.xyz*bohrrad/angstrom;  % bohr -> Angstrom

% g tensor: symmetrize raw g matrix and diagonalize
if ~isempty(d.g_matrix)
  gsym = sqrtm(d.g_matrix*d.g_matrix.');
  gsym = (gsym+gsym.')/2;  % eliminate numerical errors
  [Sys.g,Sys.gFrame] = diagonalizetensor(gsym);
end

% D tensor: convert to MHz and diagonalize
if ~isempty(d.d_raw)
  Draw = d.d_raw*100*clight/1e6;  % cm^-1 -> MHz
  Draw = (Draw+Draw.')/2;  % eliminate numerical asymmetry
  [Sys.D,Sys.DFrame] = diagonalizetensor(Draw);
end

% Hyperfine and quadrupole tensors
nAtoms = numel(d.Elements);
A = cell(1,nAtoms);
AFrame = cell(1,nAtoms);
Q = cell(1,nAtoms);
QFrame = cell(1,nAtoms);

if ~isempty(d.ATensor)
  for n = 1:numel(d.ATensor.NucIdx)
    iAtom = d.ATensor.NucIdx(n);
    el = d.Elements{iAtom};
    Araw = d.ATensor.ARaw(:,:,n);
    % Rescale from ORCA's isotope to EasySpin's hyperfine reference isotope
    gref = referenceisotope(el);
    isotope = sprintf('%d%s',round(d.ATensor.Isotope(n)),el);
    if ~isempty(gref) && d.ATensor.I(n)>0 && ~strcmp(isotope,gref.symbol)
      Araw = Araw*gref.gn/nucgval(isotope);
    end
    Asym = (Araw+Araw.')/2;
    [A{iAtom},AFrame{iAtom}] = diagonalizetensor(Asym);
  end
end

if ~isempty(d.EFGTensor)
  for n = 1:numel(d.EFGTensor.NucIdx)
    iAtom = d.EFGTensor.NucIdx(n);
    [~,qrefEl] = referenceisotope(d.Elements{iAtom});
    if isempty(qrefEl), continue; end  % element has no isotopes with I>=1

    % Get EFG tensor, diagonalize, sort by eigenvalue magnitude
    efg_SI = d.EFGTensor.V(:,:,n)*hartree/echarge/bohrrad^2; % atomic unit -> SI unit
    efg_SI = (efg_SI+efg_SI.')/2;
    [eq,QFrame{iAtom}] = diagonalizetensor(efg_SI,true);

    % Calculate quadrupole parameters and quadrupole tensor
    Qmom = qrefEl.qm*barn;
    I = qrefEl.I;
    e2qQh = echarge*Qmom*eq(3)/planck/1e6;
    K = e2qQh/(4*I)/(2*I-1);
    eta = (eq(1)-eq(2))/eq(3);
    Q{iAtom} = K*[-(1-eta),-(1+eta),2];
  end
end

% Compile nuclear data, in order of atom index
hasA = ~cellfun(@isempty,A);
hasQ = ~cellfun(@isempty,Q);
NucsIdx = find(hasA | hasQ);
if ~isempty(NucsIdx)
  nNucs = numel(NucsIdx);
  Sys.Nucs = nuclist2string(d.Elements(NucsIdx));
  Sys.NucsIdx = NucsIdx;
  if any(hasA)
    Sys.A = zeros(nNucs,3);
    Sys.AFrame = zeros(nNucs,3);
  end
  if any(hasQ)
    Sys.Q = zeros(nNucs,3);
    Sys.QFrame = zeros(nNucs,3);
  end
  for n = 1:nNucs
    iAtom = NucsIdx(n);
    if hasA(iAtom)
      Sys.A(n,:) = A{iAtom};
      Sys.AFrame(n,:) = AFrame{iAtom};
    end
    if hasQ(iAtom)
      Sys.Q(n,:) = Q{iAtom};
      Sys.QFrame(n,:) = QFrame{iAtom};
    end
  end
end

end

%-------------------------------------------------------------------------------
function [vals,angles] = diagonalizetensor(T,sortByMagnitude)
if nargin<2
  sortByMagnitude = false;
end
[R_T2M,T_] = eig(T);
vals = diag(T_).';
if sortByMagnitude
  [~,idx] = sort(abs(vals));
  vals = vals(idx);
  R_T2M = R_T2M(:,idx);
end
% Enforce right-handed frame
if det(R_T2M)<0
  R_T2M(:,1) = -R_T2M(:,1);
end
angles = eulang(R_T2M.');
end

%-------------------------------------------------------------------------------
% Parse the lines of a section (without title and $End) into a structure.
% Each entry &Name becomes a field containing a cell array with one element
% per occurrence of &Name in the section.
function P = parsesection(L)
P = struct;
isEntry = ~cellfun(@isempty,regexp(L,'^\s*&','once'));
entryLines = find(isEntry);
entryLines(end+1) = numel(L)+1;
for e = 1:numel(entryLines)-1
  k = entryLines(e);
  tok = regexp(L{k},'^\s*&(\w+)\s*(\[[^\]]*\])?\s*(.*)$','tokens','once');
  name = tok{1};
  attributes = tok{2};
  valuestr = strtrim(tok{3});
  block = L(k+1:entryLines(e+1)-1);

  type = regexp(attributes,'&Type\s*"(\w+)"','tokens','once');
  if isempty(type), type = ''; else, type = type{1}; end
  dims = str2double(regexp(attributes,'&Dim\s*\((\d+)\s*,\s*(\d+)\)','tokens','once'));

  switch type
    case {'Integer','Double',''}
      value = sscanf(valuestr,'%f',1);
    case 'Boolean'
      value = strncmpi(valuestr,'true',4);
    case 'String'
      value = regexp(valuestr,'^"([^"]*)"','tokens','once');
      value = value{1};
    case {'ArrayOfDoubles','ArrayOfIntegers'}
      value = parsearray(block,dims,name);
    case 'Coordinates'
      value = parsecoordinates(block,dims);
    otherwise
      value = valuestr;  % unknown type: keep as string
  end

  if ~isfield(P,name)
    P.(name) = {};
  end
  P.(name){end+1} = value;
end
end

%-------------------------------------------------------------------------------
% Parse an array given as blocks of columns: a header line with column
% indices, followed by lines with row index and values.
function M = parsearray(block,dims,name)
M = nan(dims);
cols = [];
for k = 1:numel(block)
  line = block{k};
  if isempty(strtrim(line)), continue; end
  if isstrprop(line(1),'digit')  % data row
    v = sscanf(line,'%f').';
    if isempty(cols) || numel(v)~=numel(cols)+1
      error('Could not parse array &%s.',name);
    end
    M(v(1)+1,cols) = v(2:end);
  else  % column header
    cols = sscanf(line,'%d').' + 1;
  end
end
if any(isnan(M(:)))
  error('Array &%s is incomplete.',name);
end
end

%-------------------------------------------------------------------------------
function value = parsecoordinates(block,dims)
nAtoms = dims(1);
block = block(~cellfun(@(x)isempty(strtrim(x)),block));
if numel(block)<nAtoms
  error('Incomplete list of Cartesian coordinates.');
end
Elements = cell(1,nAtoms);
xyz = zeros(nAtoms,3);
for iAtom = 1:nAtoms
  tok = regexp(block{iAtom},'^\s*([A-Za-z]+)\S*\s+(\S+)\s+(\S+)\s+(\S+)','tokens','once');
  if isempty(tok)
    error('Could not parse Cartesian coordinates for atom %d.',iAtom);
  end
  Elements{iAtom} = tok{1};
  xyz(iAtom,:) = str2double(tok(2:4));
end
value = {Elements,xyz};
end

%-------------------------------------------------------------------------------
% Get value of single entry &name
function varargout = getvalue(P,name,section)
if ~isfield(P,name)
  error('Entry &%s missing in section $%s.',name,section);
end
value = P.(name){1};
if nargout>1
  varargout = value;
else
  varargout = {value};
end
end

%-------------------------------------------------------------------------------
% Get values of repeated scalar entry &name as a row vector
function values = getlist(P,name,section)
if ~isfield(P,name)
  error('Entry &%s missing in section $%s.',name,section);
end
values = [P.(name){:}];
end

%-------------------------------------------------------------------------------
function checknuclearlists(T,fields,matrixfield,section)
nNucs = numel(T.NucIdx);
for f = 1:numel(fields)
  if numel(T.(fields{f}))~=nNucs
    error('Inconsistent number of nuclei in section $%s.',section);
  end
end
if size(T.(matrixfield),3)~=nNucs
  error('Inconsistent number of tensors in section $%s.',section);
end
end
