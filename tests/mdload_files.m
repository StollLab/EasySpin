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
OutOpt.keepProtCA = 1;

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
    elseif ~areequal(Traj.ProtCAxyz, ProtCAxyz_ref)...
           ||~areequal(Traj.FrameTraj, FrameTraj_ref)...
           ||~areequal(Traj.FrameTrajwrtProt, FrameTrajwrtProt_ref)...
           ||~areequal(Traj.dihedrals, dihedrals_ref)...
           ||~areequal(Traj.RProtDiff, RProtDiff_ref)
      readerr(iFile) = true;
      fprintf('   Loaded trajectories did not match reference trajectories for:\n   "%s" and "%s".\n',...
              Files{iFile,1},Files{iFile,2})
    end
  end

else

  AtomInfo.TopFile = Files{1,2};
  Traj = mdload(Files{1,1}, AtomInfo, OutOpt);
  
  ProtCAxyz_ref = Traj.ProtCAxyz;
  FrameTraj_ref = Traj.FrameTraj;
  FrameTrajwrtProt_ref = Traj.FrameTrajwrtProt;
  dihedrals_ref = Traj.dihedrals;
  RProtDiff_ref = Traj.RProtDiff;
  
  save(OldDataFile, 'ProtCAxyz_ref', 'FrameTraj_ref', 'FrameTrajwrtProt_ref',...
                    'dihedrals_ref', 'RProtDiff_ref')
  
  readerr = 0;
  
end

err = any(readerr);

data = [];

end
