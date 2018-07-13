%  mdsave  Save processed molecular dynamics simulation data to a .mat file.
%
%   mdsave(FileName, MD);
%
%   Input:
%     FileName       string
%                    Name of output file.
%
%     MD             structure array containing the following fields:
%
%                    nSteps      integer
%                                total number of steps in trajectory
%
%                    dt          double
%                                size of time step (in s)
%
%                    FrameTraj   numeric array, size = (3,3,nTraj,nSteps)
%                                xyz coordinates of coordinate frame axis
%                                vectors, x-axis corresponds to
%                                FrameTraj(:,1,nTraj,:), y-axis corresponds to
%                                FrameTraj(:,2,nTraj,:), etc.
%
%                    RProtDiff   numeric array, size = (3,3,nTraj,nSteps)
%                                trajectories of protein global rotational
%                                diffusion represented by rotation matrices
%
%             FrameTrajwrtProt   numeric array, size = (3,3,nTraj,nSteps)
%                                same as FrameTraj, but with global
%                                rotational diffusion of protein removed
%
%                    dihedrals   numeric array, size = (5,nTraj,nSteps)
%                                dihedral angles of spin label side chain
%                                bonds


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

if ~isfield(Opt, 'overwrite')
  % file overwriting is not allowed by default
  Opt.overwrite = 0;
end

if exist(FileName, 'file')
  if Opt.overwrite==0
    error('The file %s already exists. Choose a different filename or overwrite the file using Opt.overwrite.', FileName)
  else
    delete(FileName)
  end
end

if ~isfield(MD,'nSteps')
  error('The number of steps MD.nSteps must be provided.')
end

if ~isfield(MD,'dt')
  error('The time step MD.dt must be provided.')
end

if ~isfield(MD,'FrameTraj')
  error('MD.FrameTraj must be provided.')
end

if ~isfield(MD,'FrameTrajwrtProt')
  error('MD.FrameTrajwrtProt must be provided.')
end

if ~isfield(MD,'dihedrals')
  error('MD.dihedrals must be provided.')
end
  
save(FileName,MDname)

end