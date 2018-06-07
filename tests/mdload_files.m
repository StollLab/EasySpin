function [err,data] = test(opt,olddata)

BaseDir = 'mdfiles/';
SimpleFile = [BaseDir, 'A10R1_polyAla'];

load([BaseDir, 'ONxyz_ref.mat'])
load([BaseDir, 'NNxyz_ref.mat'])
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

AtomInfo.AtomNames.ONname = 'ON';
AtomInfo.AtomNames.NNname = 'NN';
AtomInfo.AtomNames.C1name = 'C1';
AtomInfo.AtomNames.C2name = 'C2';
AtomInfo.AtomNames.C1Rname = 'C1R';
AtomInfo.AtomNames.C2Rname = 'C2R';
AtomInfo.AtomNames.C1Lname = 'C1L';
AtomInfo.AtomNames.S1Lname = 'S1L';
AtomInfo.AtomNames.SGname = 'SG';
AtomInfo.AtomNames.CBname = 'CB';
AtomInfo.AtomNames.CAname = 'CA';
AtomInfo.AtomNames.Nname = 'N';

OutOpt = [];

nTests = size(Files, 1);

for iFile = 1:nTests
  readerr(iFile) = false;
  AtomInfo.TopFile = Files{iFile,2};
  Traj = mdload(Files{iFile,1}, AtomInfo, OutOpt);
  if any(~structfun(@(x) areequal(isnan(x),0), Traj))
    readerr(iFile) = true;
    fprintf('   NaNs were detected in output from mdload with args:\n   "%s" and "%s".\n',...
            Files{iFile,1},Files{iFile,2})
  elseif ~areequal(Traj.ONxyz, ONxyz_ref)||~areequal(Traj.NNxyz, NNxyz_ref)...
         ||~areequal(Traj.C1xyz, C1xyz_ref)||~areequal(Traj.C2xyz, C2xyz_ref)
    readerr(iFile) = true;
    fprintf('   Loaded trajectories did not match reference trajectories for:\n   "%s" and "%s".\n',...
            Files{iFile,1},Files{iFile,2})
  end
end

err = any(readerr);

data = [];

end
