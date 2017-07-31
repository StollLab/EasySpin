function [err,data] = test(opt,olddata)

BaseDir = 'mdfiles/';
SimpleFile = [BaseDir, 'A10R1_polyAla'];

load([BaseDir, 'Oxyz_ref.mat'])
load([BaseDir, 'Nxyz_ref.mat'])
load([BaseDir, 'C1xyz_ref.mat'])
load([BaseDir, 'C2xyz_ref.mat'])

%-------------------------------------------------
% Read several formats
%-------------------------------------------------

Extensions = {...
'.dcd','.psf';'.DCD','.PSF'
};

Files = strcat(SimpleFile, Extensions);

AtomInfo.ResName = 'CYR1';

AtomInfo.AtomNames.OName = 'ON';
AtomInfo.AtomNames.NName = 'NN';
AtomInfo.AtomNames.C1Name = 'C1';
AtomInfo.AtomNames.C2Name = 'C2';

OutOpt.Frame = 0;

nTests = size(Files, 1);

for iFile = 1:nTests
  readerr(iFile) = false;
  AtomInfo.TopFile = Files{iFile,2};
  Traj = mdload(Files{iFile,1}, AtomInfo, OutOpt);
  if any(~structfun(@(x) areequal(isnan(x),0), Traj))
    readerr(iFile) = true;
    fprintf('   NaNs were detected in output from mdload with args:\n   "%s" and "%s".\n',...
            Files{iFile,1},Files{iFile,2})
  elseif ~areequal(Traj.Oxyz, Oxyz_ref)||~areequal(Traj.Nxyz, Nxyz_ref)...
         ||~areequal(Traj.C1xyz, C1xyz_ref)||~areequal(Traj.C2xyz, C2xyz_ref)
    readerr(iFile) = true;
    fprintf('   Loaded trajectories did not match reference trajectories for:\n   "%s" and "%s".\n',...
            Files{iFile,1},Files{iFile,2})
  end
end

err = any(readerr);

data = [];

end
