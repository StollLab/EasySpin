function [err,data] = test(opt,olddata)

BaseDir = 'mdfiles/';
SimpleFile = [BaseDir, 'A10R1_polyAla'];

%-------------------------------------------------
% Read several formats
%-------------------------------------------------

Extensions = {...
'.dcd','.psf';'.DCD','.PSF'
};

Files = strcat(SimpleFile, Extensions);

AtomInfo.ResName = 'CYR1';

AtomInfo.SegName = 'PROT';

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

OutOpt.Verbosity = 0;

nTests = size(Files, 1);

OldDataFile = [BaseDir, 'mdload_A10R1.mat'];

if exist(OldDataFile,'file')>0
  
  load(OldDataFile)

  for iFile = 1:nTests
    readerr(iFile) = false;
    AtomInfo.TopFile = Files{iFile,2};
    Traj = mdload(Files{iFile,1}, AtomInfo, OutOpt);
    if any(~structfun(@(x) areequal(isnan(x),0), Traj))
      readerr(iFile) = true;
      fprintf('   NaNs were detected in output from mdload with args:\n   "%s" and "%s".\n',...
              Files{iFile,1},Files{iFile,2})
    elseif ~areequal(Traj.Protxyz, Protxyz_ref)||~areequal(Traj.Labelxyz, Labelxyz_ref)...
           ||~areequal(Traj.FrameZ, FrameZ_ref)||~areequal(Traj.chi1, chi1_ref)
      readerr(iFile) = true;
      fprintf('   Loaded trajectories did not match reference trajectories for:\n   "%s" and "%s".\n',...
              Files{iFile,1},Files{iFile,2})
    end
  end

else

  AtomInfo.TopFile = Files{1,2};
  Traj = mdload(Files{1,1}, AtomInfo, OutOpt);
  
  Protxyz_ref = Traj.Protxyz;
  Labelxyz_ref = Traj.Labelxyz;
  FrameZ_ref = Traj.FrameZ;
  chi1_ref = Traj.chi1;
  
  save(OldDataFile, 'Protxyz_ref', 'Labelxyz_ref', 'FrameZ_ref', 'chi1_ref')
  
  readerr = 0;
  
end

err = any(readerr);

data = [];

end
