%  readdcd  Read a binary DCD file obtained from a molecular dynamics
%           simulation.
%
%   traj = readdcd(fileName);
%   traj = readdcd(fileName, idx);
%
%   Input:
%     fileName       character array
%                    Name of DCD file to be read.
%
%     idx            numeric
%                    Index or indices of atoms to be read.
%
%   Output:
%     traj           numeric, size = (nFrames,3)
%                    x,y,z coordinates of desired atoms in the trajectory

% current code is based on 'readdcd' in MDToolbox:
%  https://github.com/ymatsunaga/mdtoolbox

function traj = readdcd(varargin)

switch nargin
  case 1
    fileName = varargin{1};
    idx = [];
  case 2
    fileName = varargin{1};
    idx = varargin{2};
    nidx = length(idx);
end

if ~ischar(fileName)
  error('Filename must be given as a character array.')
end

if ~endsWith(fileName,'.dcd')
  error('Please give the full filename including the ".dcd" extension.')
end

if islogical(idx)
  idx = find(idx);
end

% obtain integer file identifier
fileID = fopen(fileName, 'r');

% check size of file
fseek(fileID, 0, 'eof');
fileSize = ftell(fileID);
fseek(fileID, 0, 'bof');

% check endianness of binary file
% TODO can we just check machineformat and be done?
[fileName, ~, machineformat] = fopen(fileID);

% block 1
% -------------------------------------------------------------------------
header.blocksize1 = fread(fileID, 1, 'int32');

% TODO figure out why this needs to be 84
if header.blocksize1 ~= 84
  fclose(fileID);
  error('The following file does not have proper format:\n"%s"', fileName)
end

% CORD for coordinates or VELD for velocities
header.HDR = fread(fileID, 4, 'char');

% number of frames
header.NSET = fread(fileID, 1, 'int32');

% starting time step
header.ISTRT = fread(fileID, 1, 'int32');

% trajectory saving frequency
header.NSAVC = fread(fileID, 1, 'int32');

% total number of time steps
header.NSTEP = fread(fileID, 1, 'int32');

% null4 (int*4)
header.NULL4 = fread(fileID, 4, 'int32');

% number of free atoms
header.NFREAT = fread(fileID, 1, 'int32');

% size of time step
header.DELTA = fread(fileID, 1, 'float32');

% null9 (int*9)
header.NULL9 = fread(fileID, 9, 'int32');

% version
header.VERSION = fread(fileID, 1, 'int32');

% check for charmm format
if header.VERSION > 0
  header.isCHARMM = 1;
  
  cof = ftell(fileID);

  % CHARMM extrablock
  fseek(fileID, 48, 'bof');
  if fread(fileID, 1, 'int32')==1
    header.CHARMMextrablock = 1;
  end
  
  % CHARMM 4dims
  fseek(fileID, 52, 'bof');
  if fread(fileID, 1, 'int32')==1
    header.CHARMMextrablock = 1;
  end

  fseek(fileID, cof, 'bof');

else
  % xplor format
  header.isCHARMM = 0;
  
  cof = ftell(fileID);
  fseek(fileID, 44, 'bof');
  header.DELTA = fread(fileID, 1, 'float64');
  fseek(fileID, cof, 'bof');
end

% check consistency
blocksize1 = fread(fileID, 1, 'int32');
if header.blocksize1~=blocksize1
  fclose(fileID);
  error('Inconsistent sizes for block 1.')
end

% block 2 (title)
% -------------------------------------------------------------------------
% blocksize2
header.blocksize2 = fread(fileID, 1, 'int32');

% # of title lines
header.ntitle = fread(fileID, 1, 'int32');

% title
header.TITLE = fread(fileID, 80*header.ntitle, 'char');
header.TITLE = char(reshape(header.title, 80, [])');

% check consistency
blocksize2 = fread(fileID, 1, 'int32');
if header.blocksize2~=blocksize2
  fclose(fileID);
  error('Inconsistent sizes for block 2.')
end

% block 3 (natom)
% -------------------------------------------------------------------------
% blocksize3
header.blocksize3 = fread(fileID, 1, 'int32');

% # of atoms
header.NATOM = fread(fileID, 1, 'int32');

% check consistency
blocksize3 = fread(fileID, 1, 'int32');
if header.blocksize3~=blocksize3
  fclose(fileID);
  error('Inconsistent sizes for block 3.')
end

% coordinates
% -------------------------------------------------------------------------
headersize = ftell(fileID);

if header.CHARMMextrablock
  extrablocksize = 4*2 + 8*6;
else
  extrablocksize = 0;
end

coordblocksize = (4*2 + 4*header.NATOM)*3;

nFrames = floor(fileSize - headersize) / (extrablocksize + coordblocksize);

if isempty(idx)
  idx = 1:header.NATOM;
end

traj = zeros(nFrames, numel(idx)*3);
box = zeros(nFrames, 3);

% read next frames
for iFrame = 1:nFrames
  % charmm extrablock (unitcell info)
  if header.CHARMMextrablock
    blocksize = fread(fileID, 1, 'int32');
    dummy = fread(fileID, 6, 'float64');
    blocksize = fread(fileID, 1, 'int32');
  end

  % x coordinates
  blocksize = fread(fileID, 1, 'int32');
  x = fread(fileID, blocksize/4, 'float32');
  blocksize = fread(fileID, 1, 'int32');

  % y coordinates 
  blocksize = fread(fileID, 1, 'int32');
  y = fread(fileID, blocksize/4, 'float32');
  blocksize = fread(fileID, 1, 'int32');

  % z coordinates 
  blocksize = fread(fileID, 1, 'int32');
  z = fread(fileID, blocksize/4, 'float32');
  blocksize = fread(fileID, 1, 'int32');
  
%   % ignore charmm 4dims extension
%   if header.CHARMM4dims
%     blocksize = fread(fileID, 1, 'int32');
%     fseek(fileID, blocksize, 0);
%     blocksize = fread(fileID, 1, 'int32');
%   end

%   if header.CHARMMextrablock
%     box(iFrame, :) = dummy([1 3 6])';
%   end
  traj(iFrame, 1:3:end) = x(idx)';
  traj(iFrame, 2:3:end) = y(idx)';
  traj(iFrame, 3:3:end) = z(idx)';
end
  
fclose(fileID);

end