function [err,data] = test(opt,olddata)

BaseDir = 'mdfiles/';
SimpleFile = [BaseDir, 'A10TOAC_polyAla'];

%-------------------------------------------------
% Read several formats
%-------------------------------------------------

Extensions = {...
'.dcd','.psf';'.DCD','.PSF'
};

Files = strcat(SimpleFile, Extensions);

AtomInfo.ResName = 'TOC';

AtomInfo.SegName = 'PROA';

AtomInfo.LabelName = 'TOAC';

AtomInfo.AtomNames.ONname = 'OE';
AtomInfo.AtomNames.NNname = 'ND';
AtomInfo.AtomNames.CGSname = 'CG2';
AtomInfo.AtomNames.CGRname = 'CG1';
AtomInfo.AtomNames.CBSname = 'CB2';
AtomInfo.AtomNames.CBRname = 'CB1';
AtomInfo.AtomNames.CAname = 'CA';
AtomInfo.AtomNames.Nname = 'N';

OutOpt.Verbosity = 0;
OutOpt.keepProtCA = 1;

OutOpt.Verbosity = 0;
OutOpt.keepProtCA = 1;

AtomInfo.TopFile = Files{1,2};
data = mdload(Files{1,1}, AtomInfo, OutOpt);

% Compare
% -------------------------------------------------------------------------

nTests = size(Files, 1);

if ~isempty(olddata)
  for iFile = 1:nTests
    err(iFile) = false;
    AtomInfo.TopFile = Files{iFile,2};
    data = mdload(Files{iFile,1}, AtomInfo, OutOpt);
%     if any(~structfun(@(x) areequal(isnan(x),0), Traj))
%       readerr(iFile) = true;
%       fprintf('   NaNs were detected in output from mdload with args:\n   "%s" and "%s".\n',...
%               Files{iFile,1},Files{iFile,2})
    if ~areequal(data.ProtCAxyz, olddata.ProtCAxyz)...
           ||~areequal(data.FrameTraj, olddata.FrameTraj)...
           ||~areequal(data.FrameTrajwrtProt, olddata.FrameTrajwrtProt)...
           ||~areequal(data.dihedrals, olddata.dihedrals)...
           ||~areequal(data.RProtDiff, olddata.RProtDiff)
      err(iFile) = true;
      fprintf('   Loaded trajectories did not match reference trajectories for:\n   "%s" and "%s".\n',...
              Files{iFile,1},Files{iFile,2})
    end
  end
  err = any(err);
else
  err = [];  
end

end
