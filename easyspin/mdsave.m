%  mdsave  Save processed molecular dynamics simulation data.
%
%   mdsave(FileName, Traj, Opt);
%
%   Input:
%     FileName       character array
%                    Name of output file.
%
%     MD             structure array containing the following fields:
%
%                    nSteps   integer
%                             total number of steps in trajectory
%
%                    dt       double
%                             size of time step (in s)
%
%                    if Frame = 0:
%                      ONxyz     numeric, size = (nSteps,3)
%                               x,y,z coordinate trajectory for oxygen
% 
%                      NNxyz     numeric, size = (nSteps,3)
%                               x,y,z coordinate trajectory for nitrogen
%  
%                      C1xyz    numeric, size = (nSteps,3)
%                               x,y,z coordinate trajectory for carbon C1
%  
%                      C2xyz    numeric, size = (nSteps,3)
%                               x,y,z coordinate trajectory for carbon C2
%
%                    if Frame = 1;
%                      FrameX   numeric, size = (nSteps,3)
%                               x,y,z coordinate trajectory for X-axis
%                               vector
% 
%                      FrameY   numeric, size = (nSteps,3)
%                               x,y,z coordinate trajectory for Y-axis
%                               vector
%  
%                      FrameZ   numeric, size = (nSteps,3)
%                               x,y,z coordinate trajectory for Z-axis
%                               vector
%
%     Opt            structure array containing the following fields:
%
%
%   Supported formats are identified via the extension
%   in 'TrajFile' and 'TopFile'. Extensions:
%
%     NAMD, CHARMM:        .DCD, .PSF
%

function mdsave(FileName, MD, Opt)

switch nargin
  case 0
    help(mfilename); return;
  case 2 % FileName and MD specified, initialize Opt
    Opt = struct;
  case 3 % all inputs given
  otherwise
    error('Incorrect number of input arguments.')
end

% store name of variable given for MD argument
MDname = inputname(2);

if ~isfield(Opt,'Frame'), Opt.Frame = 1; end
if ~isfield(Opt, 'overwrite')
  % file overwriting is not allowed by default
  Opt.overwrite = 0;
end

if exist(FileName, 'file')
  if Opt.overwrite==0
    error('%s already exists. Choose a different filename or enable overwriting.', FileName)
  else
    delete(FileName)
  end
end

% % supported file types
% supportedTrajFileExts = {'.DCD'};
% supportedTopFileExts = {'.PSF'};
% 
% TopFile = AtomInfo.TopFile;
% ResName = AtomInfo.ResName;
% AtomNames = AtomInfo.AtomNames;
% 
% if ~ischar(TopFile)||regexp(TopFile,'\w+\.\w+','once')<1
%   error('TopFile must be given as a character array, including the filename extension.')
% end
% 
% if numel(regexp(TopFile,'\.'))>1
%   error('Only one period (".") can be included in TopFile as part of the filename extension. Remove the others.')
% end
% 
% [TopFilePath, TopFileName, TopFileExt] = fileparts(TopFile);
% TopFile = fullfile(TopFilePath, [TopFileName, TopFileExt]);

if ~isfield(MD,'dt')
  error('The time step dt must be given in MD.')
end

if Opt.Frame==1
  % give the reference frame coordinate axis vector trajectories as output
  
  MD.isFrame = 1;
  
  if ~isfield(MD,'FrameTraj')
    error('If saving a frame trajectory, input MD must contain FrameTraj coordinates.')
  end
  
  if any(any(size(MD.FrameTraj)~=[3,3,1,MD.nSteps]))
    error('MD.FrameTraj must be of size (3,3,1,MD.nSteps).')
  end
  
end
  
save(FileName,MDname)

end