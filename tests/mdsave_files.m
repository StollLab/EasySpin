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


OutOpt.Type = 'Frame';
OutOpt.overwrite = 1;

nTests = size(Files, 1);

tempFileName = [BaseDir, 'temp'];

for iFile = 1:nTests
  readerr(iFile) = false;
  AtomInfo.TopFile = Files{iFile,2};
  MD = mdload(Files{iFile,1}, AtomInfo, OutOpt);
  
  FrameX_ref = MD.FrameX;
  FrameY_ref = MD.FrameY;
  FrameZ_ref = MD.FrameZ;
  
  mdsave(tempFileName,MD,OutOpt)
  load(tempFileName)
  
  if ~areequal(MD.FrameX, FrameX_ref)||~areequal(MD.FrameY, FrameY_ref)...
         ||~areequal(MD.FrameZ, FrameZ_ref)
    readerr(iFile) = true;
    fprintf('   Saved trajectories did not match reference trajectories for:\n   "%s" and "%s".\n',...
            Files{iFile,1},Files{iFile,2})
  end
end

err = any(readerr);

data = [];

delete([tempFileName, '.mat'])

end
