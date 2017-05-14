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

OutOpt.Frame = 1;
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

err = readerr;

data = [];

delete([tempFileName, '.mat'])

end
