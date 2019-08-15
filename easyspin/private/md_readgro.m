% md_readgro  Read protein structure file in Gromos87 format.
%
%   data = md_readgro(FileName);
%
% Input:
%   grofile    filename, including .gro extension
%
% Output:
%   data      structure with field containing all data read from the file


function data = md_readgro(GroFile,ResName,LabelName,AtomNames)

if ~strcmpi(GroFile(end-3:end),'.GRO')
  error('File extension must be .gro.');
end

fh = fopen(GroFile);
if fh<0
  error('Could not open file ''%s''',GroFile);
end

% Read entire file
allLines = textscan(fh,'%s','whitespace','','delimiter','\n');
fclose(fh);
allLines = allLines{1};

title = allLines{1};
nAtoms = sscanf(allLines{2},'%d');

if numel(allLines)~=nAtoms+3
  error('Number of lines in file ''%s'' appears incorrect.',GroFile);
end

resnum = zeros(1,nAtoms);
resnames = cell(1,nAtoms);
atomnames = cell(1,nAtoms);
atomnumber = zeros(1,nAtoms);
pos = zeros(nAtoms,3);
vel = zeros(nAtoms,3);
for k = 1:nAtoms
  L = allLines{k+2};
  resnum(k) = sscanf(L(1:5),'%d');
  resnames{k} = L(6:10);
  atomnames{k} = L(11:15);
  atomnumber(k) = sscanf(L(16:20),'%d');
  pos(k,:) = sscanf(L(21:44),'%f %f %f');
  if length(L)>44
    vel(k,:) = sscanf(L(45:68),'%f %f %f');
  end
end
box = sscanf(allLines{nAtoms+3},'%f');
resnames = strtrim(resnames);
atomnames = strtrim(atomnames);

data.idx_ProteinCA = strcmpi(atomnames,'CA');
idx_SpinLabel = strcmpi(resnames,ResName);
data.idx_SpinLabel = idx_SpinLabel;

% Locate spin label atoms
findatomindex = @(name) strcmpi(atomnames(idx_SpinLabel),name);
switch LabelName
  case 'R1'
    data.idx_ON = findatomindex(AtomNames.ONname);
    data.idx_NN = findatomindex(AtomNames.NNname);
    data.idx_C1 = findatomindex(AtomNames.C1name);
    data.idx_C2 = findatomindex(AtomNames.C2name);
    data.idx_C1R = findatomindex(AtomNames.C1Rname);
    data.idx_C2R = findatomindex(AtomNames.C2Rname);
    data.idx_C1L = findatomindex(AtomNames.C1Lname);
    data.idx_S1L = findatomindex(AtomNames.S1Lname);
    data.idx_SG = findatomindex(AtomNames.SGname);
    data.idx_CB = findatomindex(AtomNames.CBname);
    data.idx_CA = findatomindex(AtomNames.CAname);
    data.idx_N = findatomindex(AtomNames.Nname);
  case 'TOAC'
    data.idx_ON = findatomindex(AtomNames.ONname);
    data.idx_NN = findatomindex(AtomNames.NNname);
    data.idx_CGS = findatomindex(AtomNames.CGSname);
    data.idx_CGR = findatomindex(AtomNames.CGRname);
    data.idx_CBS = findatomindex(AtomNames.CBSname);
    data.idx_CBR = findatomindex(AtomNames.CBRname);
    data.idx_CA = findatomindex(AtomNames.CAname);
    data.idx_N = findatomindex(AtomNames.Nname);
end

% Assemble output structure
data.title = title;
data.nAtoms = nAtoms;
data.resnum = resnum;
data.resname = resnames;
data.atomname = atomnames;
data.atomnumber = atomnumber;
data.pos = pos;
data.vel = vel;
data.box = box;

return
