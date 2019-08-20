%  md_readdcd  Read a binary DCD file obtained from a molecular dynamics
%                    simulation.
%
%   traj = md_readdcd(FileName);
%   traj = md_readdcd(FileName, idx);
%
%   Input:
%     FileName       character array
%                    Name of DCD file to be read.
%
%     idx            numeric
%                    Index or indices of atoms to be read.
%
%
%   Output:
%     Traj           structure with trajectory data
%        .nAtoms     number of atoms
%        .nFrames    total number of frames
%        .xyz        x,y,z coordinates, 3D-array, (nFrames,3,nAtoms)
%        .dt         time step, in seconds

% current code is based on 'readdcd' in MDToolbox:
%  https://github.com/ymatsunaga/mdtoolbox

function Traj = md_readdcd(varargin)

switch nargin
  case 1
    FileName = varargin{1};
    idx = [];
  case 2
    FileName = varargin{1};
    idx = varargin{2};
end

if ~ischar(FileName)
  error('Filename must be given as a character array.')
end

if ~strcmpi(FileName(end-3:end),'.dcd')
  error('Please give the full filename including the ".dcd" extension.')
end

if islogical(idx)
  idx = find(idx);
end

% obtain integer file identifier
FileID = fopen(FileName, 'r');
if FileID<1, error('File "%s" could not be opened.', FileName); end
ensurefclose = onCleanup(@() fclose(FileID));

% check size of file
fseek(FileID, 0, 'eof');
fileSize = ftell(FileID);
fseek(FileID, 0, 'bof');

% check endianness of binary file
[FileName, ~, machineformat] = fopen(FileID);

% block 1
% -------------------------------------------------------------------------
header.blocksize1 = fread(FileID, 1, 'int32');

% TODO figure out why this needs to be 84
if header.blocksize1 ~= 84
  error('File "%s" does not have proper format.', FileName)
end

% CORD for coordinates or VELD for velocities
header.HDR = fread(FileID, 4, 'char');

% number of frames
header.NSET = fread(FileID, 1, 'int32');

% starting time step
header.ISTRT = fread(FileID, 1, 'int32');

% trajectory saving frequency
header.NSAVC = fread(FileID, 1, 'int32');

% total number of time steps
header.NSTEP = fread(FileID, 1, 'int32');

% null4 (int*4)
header.NULL4 = fread(FileID, 4, 'int32');

% number of free atoms
header.NFREAT = fread(FileID, 1, 'int32');

% size of time step
header.DELTA = fread(FileID, 1, 'float32');

% null9 (int*9)
header.NULL9 = fread(FileID, 9, 'int32');

% version
header.VERSION = fread(FileID, 1, 'int32');

% check for charmm format
if header.VERSION > 0
  header.isCHARMM = 1;
  
  cof = ftell(FileID);

  % CHARMM extrablock
  fseek(FileID, 48, 'bof');
  if fread(FileID, 1, 'int32')==1
    header.CHARMMextrablock = 1;
  end
  
  % CHARMM 4dims
  fseek(FileID, 52, 'bof');
  if fread(FileID, 1, 'int32')==1
    header.CHARMMextrablock = 1;
  end

  fseek(FileID, cof, 'bof');

else
  % xplor format
  header.isCHARMM = 0;
  
  cof = ftell(FileID);
  fseek(FileID, 44, 'bof');
  header.DELTA = fread(FileID, 1, 'float64');
  fseek(FileID, cof, 'bof');
end

% check consistency
blocksize1 = fread(FileID, 1, 'int32');
if header.blocksize1~=blocksize1
  error('Inconsistent sizes for block 1.')
end

% block 2 (title)
% -------------------------------------------------------------------------
% blocksize2
header.blocksize2 = fread(FileID, 1, 'int32');

% # of title lines
header.NTITLE = fread(FileID, 1, 'int32');

% title
header.TITLE = fread(FileID, 80*header.NTITLE, 'char');
header.TITLE = char(reshape(header.TITLE, 80, [])');

% check consistency
blocksize2 = fread(FileID, 1, 'int32');
if header.blocksize2~=blocksize2
  error('Inconsistent sizes for block 2.')
end

% block 3 (natom)
% -------------------------------------------------------------------------
% blocksize3
header.blocksize3 = fread(FileID, 1, 'int32');

% # of atoms
header.NATOM = fread(FileID, 1, 'int32');

% check consistency
blocksize3 = fread(FileID, 1, 'int32');
if header.blocksize3~=blocksize3
  error('Inconsistent sizes for block 3. Check the integrity of the DCD file.')
end

% coordinates
% -------------------------------------------------------------------------
headersize = ftell(FileID);

if header.CHARMMextrablock
  extrablocksize = 4*2 + 8*6;
else
  extrablocksize = 0;
end

coordblocksize = (4*2 + 4*header.NATOM)*3;

nSnapShots = floor(fileSize - headersize) / (extrablocksize + coordblocksize);

% if nSnapShots~=header.NSTEP-header.ISTRT  FIXME why doesn't this work?
%   error('Number of snapshots is not equal to the number of time steps. Check the integrity of the DCD file.')
% end

if isempty(idx)
  idx = 1:header.NATOM;
end

Traj.xyz = zeros(nSnapShots, numel(idx)*3);
box = zeros(nSnapShots, 3);

% read next frames
for iFrame = 1:nSnapShots
  % charmm extrablock (unitcell info)
  if header.CHARMMextrablock
    blocksize = fread(FileID, 1, 'int32');
    dummy = fread(FileID, 6, 'float64');
    blocksize = fread(FileID, 1, 'int32');
  end

  % x coordinates
  blocksize = fread(FileID, 1, 'int32');
  x = fread(FileID, blocksize/4, 'float32');
  blocksize = fread(FileID, 1, 'int32');

  % y coordinates 
  blocksize = fread(FileID, 1, 'int32');
  y = fread(FileID, blocksize/4, 'float32');
  blocksize = fread(FileID, 1, 'int32');

  % z coordinates 
  blocksize = fread(FileID, 1, 'int32');
  z = fread(FileID, blocksize/4, 'float32');
  blocksize = fread(FileID, 1, 'int32');
  
%   % ignore charmm 4dims extension
%   if header.CHARMM4dims
%     blocksize = fread(fileID, 1, 'int32');
%     fseek(fileID, blocksize, 0);
%     blocksize = fread(fileID, 1, 'int32');
%   end

%   if header.CHARMMextrablock
%     box(iFrame, :) = dummy([1 3 6])';
%   end
  Traj.xyz(iFrame, 1:3:end) = x(idx)';
  Traj.xyz(iFrame, 2:3:end) = y(idx)';
  Traj.xyz(iFrame, 3:3:end) = z(idx)';
end
Traj.xyz = reshape(Traj.xyz, [nSnapShots, 3, numel(idx)]);

isNAMD = true;
if isNAMD
  timeunit = 48.88821; % NAMD internal time unit, in femtoseconds
else
  timeunit = 1; % femtoseconds
end
Traj.dt = header.NSAVC*header.DELTA*timeunit*1e-15; % femtoseconds -> seconds

Traj.nAtoms = header.NATOM;
Traj.nFrames = nSnapShots;

end
