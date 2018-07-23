%  md_savedcd  Save a binary DCD file containing processed molecular dynamics 
%                    simulation data.
%
%   md_savedcd(FileName, Traj);
%
%   Input:
%     FileName       character array
%                    Name of DCD file to be read.
%
%     Traj           structure array
%                    xyz: x,y,z coordinates of desired atoms in the trajectory
%                    nSteps: total number of time steps in the trajectory
%                    dt: size of the time step (in sec)
%
%     Opt            structure array
%                    overwrite: 0 throw error if FileName already exists
%                               1 overwrite existing FileName
%

% current code is based on 'savedcd' in MDToolbox:
%  https://github.com/ymatsunaga/mdtoolbox

function md_savedcd(FileName, Traj, Opt)

switch nargin
  case 2
    % Opt not given
    Opt = struct();
  case 3
    % everything given
  otherwise
    error('Wrong number of input arguments.')
end

if ~ischar(FileName)
  error('Filename must be given as a character array.')
end

if ~endsWith(lower(FileName),'.dcd')
  error('Please give the full filename including the ".dcd" extension.')
end

if ~isstruct(Traj)
  error('Traj must be a structured array.')
end

if ~isfield(Opt, 'overwrite')
  % file overwriting is not allowed by default
  Opt.overwrite = 0;
end

if exist(FileName, 'file')&&Opt.overwrite==0
  error('%s already exists. Choose a different filename or enable overwriting.', FileName)
end

if Opt.overwrite==1
  delete FileName
end

% obtain integer file identifier
FileID = fopen(FileName, 'w');
if FileID<1, error('File "%s" could not be opened.', FileName); end
ensurefclose = onCleanup(@() fclose(FileID));

TIMEFACTOR = 48.88821;  % used to convert internal time units to fs TODO this should be restricted to NAMD output

% initialize
% -------------------------------------------------------------------------
[nframe, natom3] = size(trj);
natom = natom3 / 3;

% header

% block sizes
header.blocksize1 = 84;
header.blocksize2 = 164;
header.blocksize3 = 4;

% default header in xplor format
header.is_charmm = false;
header.is_charmm_extrablock = false;
header.is_charmm_4dims = false;
header.HDR = 'CORD';
header.NSET = size(trj, 1);
header.ISTRT = 0;
header.NSAVC = 1;
header.NSTEP = 0;
header.NULL4 = zeros(4, 1);
header.NFREAT = 0;
header.DELTA = 1.0;
header.NULL9 = zeros(9, 1);
header.VERSION = 0;
header.NTITLE = 2;
title1 = sprintf('REMARKS FILENAME=%s CREATED BY MATLAB', filename);
for i = (numel(title1)+1):80
  title1 = [title1 ' '];
end
title1 = title1(1:80);
title2 = sprintf('REMARKS DATE: %s CREATED BY USER: %s', datestr(now, 'mm/dd/yy'), getenv('USER'));
for i = (numel(title2)+1):80
  title2 = [title2 ' '];
end
title2 = title2(1:80);
header.title = [title1; title2];
header.natom = size(trj, 2) / 3;

if header.nset ~= size(trj, 1)
  header.nset = size(trj, 1);
end

if exist('box_', 'var') && ~isempty(box_)
  % charmm format
  header.is_charmm = true;
  header.is_charmm_extrablock = true;
  if header.version == 0
    header.version = 1; % is_charmm -> true
  end
  header.null9(1) = 1; % is_charmm_extrablock -> true
else
  % xplor format
  header.is_charmm = false;
  header.is_charmm_extrablock = false;
  header.is_charmm_4dims = false;
  header.version = 0; % is_charmm -> false
  header.null9(1) = 0; % is_charmm_extrablock -> false
  header.null9(2) = 0; % is_charmm_4dims -> false
end

% block 1
% -------------------------------------------------------------------------

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
  error('Inconsistent sizes for block 3.')
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

nFrames = floor(fileSize - headersize) / (extrablocksize + coordblocksize);

% if nFrames~=header.NSTEP-header.ISTRT  FIXME why doesn't this work?
%   error('Number of frames is not equal to the number of time steps. Check the integrity of the DCD file.')
% end

if isempty(idx)
  idx = 1:header.NATOM;
end

Traj.xyz = zeros(nFrames, numel(idx)*3);
box_ = zeros(nFrames, 3);

% read next frames
for iFrame = 1:nFrames
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
%     box_(iFrame, :) = dummy([1 3 6])';
%   end
  Traj.xyz(iFrame, 1:3:end) = x(idx)';
  Traj.xyz(iFrame, 2:3:end) = y(idx)';
  Traj.xyz(iFrame, 3:3:end) = z(idx)';
end

Traj.dt = header.NSAVC*TIMEFACTOR*header.DELTA*1e-15;
Traj.nSteps = nFrames;
Traj.xyz = reshape(Traj.xyz, [Traj.nSteps, 3, numel(idx)]);

end