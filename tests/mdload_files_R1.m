function [ok,data] = test(opt,olddata)

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

OutOpt.Verbosity = 0;
OutOpt.keepProtCA = 1;

data = mdload(Files{1,1}, Files{1,2}, AtomInfo, OutOpt);

% Compare
% -------------------------------------------------------------------------

nTests = size(Files, 1);

thr = 1e-10;

if ~isempty(olddata)
  for iFile = 1:nTests
    data = mdload(Files{iFile,1}, Files{iFile,2}, AtomInfo, OutOpt);
    ok(iFile) = areequal(data.ProtCAxyz, olddata.ProtCAxyz,thr,'abs') && ...
       areequal(data.FrameTraj, olddata.FrameTraj,thr,'abs') && ...
       areequal(data.FrameTrajwrtProt, olddata.FrameTrajwrtProt,thr,'abs') && ...
       areequal(data.dihedrals, olddata.dihedrals,thr,'abs') && ...
       areequal(data.RProtDiff, olddata.RProtDiff,thr,'abs');
    if ~ok(iFile)
      fprintf('   Loaded trajectories did not match reference trajectories for:\n   "%s" and "%s".\n',...
              Files{iFile,1},Files{iFile,2})
    end
  end
else
  ok = [];
end

end
