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

OutOpt.Verbosity = 0;
OutOpt.keepProtCA = 1;

OutOpt.Verbosity = 0;
OutOpt.keepProtCA = 1;

data = mdload(Files{1,1}, Files{1,2}, AtomInfo, OutOpt);

% Compare
% -------------------------------------------------------------------------

nTests = size(Files, 1);

thr = 1e-10;

if ~isempty(olddata)
  for iFile = 1:nTests
    err(iFile) = false;
    data = mdload(Files{iFile,1}, Files{iFile,2}, AtomInfo, OutOpt);
%     if any(~structfun(@(x) areequal(isnan(x),0), Traj))
%       readerr(iFile) = true;
%       fprintf('   NaNs were detected in output from mdload with args:\n   "%s" and "%s".\n',...
%               Files{iFile,1},Files{iFile,2})
    if ~areequal(data.ProtCAxyz, olddata.ProtCAxyz,thr,'abs') || ...
       ~areequal(data.FrameTraj, olddata.FrameTraj,thr,'abs') || ...
       ~areequal(data.FrameTrajwrtProt, olddata.FrameTrajwrtProt,thr,'abs') || ...
       ~areequal(data.dihedrals, olddata.dihedrals,thr,'abs') || ...
       ~areequal(data.RProtDiff, olddata.RProtDiff,thr,'abs')
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
