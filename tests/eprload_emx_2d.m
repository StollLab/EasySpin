function [err,data] = eprload_esp_xaxis(opt,olddata)

data = [];
err = [];

BaseDir = 'eprfiles/';

%-------------------------------------------------
% Check whether ESP abscissa vectors are correct
%-------------------------------------------------

Files{1} = 'EMX_2Dpowersweep.par';
Files{2} = 'EMX_Acquisit_power2d.par';

XX{1} = [3452.5 75]; %XXLB/XXWI
XY{1} = [23 -24]; %XYLB/XYWI
XX{2} = [3410 140]; %XXLB/XXWI
XY{2} = [44 -36]; %XYLB/XYWI

ax_err = zeros(1,numel(Files));
for iFile = 1:numel(Files)
  if (opt.Display), disp(Files{iFile}); end
  [ax,data] = eprload([BaseDir Files{iFile}]);
  XXLB = ax{1}(1);
  XXWI = ax{1}(end)-XXLB;
  XYLB = ax{2}(1);
  XYWI = ax{2}(end)-XYLB;
  if ~areequal(XX{iFile}(1),XXLB,1e-4) || ...
     ~areequal(XX{iFile}(2),XXWI,1e-4)
    disp([' ' Files{iFile}]);
    fprintf('   x read:      %0.9f %0.9f\n',XXLB,XXWI);
    fprintf('   x should be: %0.9f %0.9f\n',XX{iFile}(1),XX{iFile},2);
    ax_err(iFile) = 1;
  end
  if ~areequal(XY{iFile}(1),XYLB,1e-4) || ...
     ~areequal(XY{iFile}(2),XYWI,1e-4)
    disp([' ' Files{iFile}]);
    fprintf('   x read:      %0.9f %0.9f\n',XYLB,XYWI);
    fprintf('   x should be: %0.9f %0.9f\n',XY{iFile}(1),XY{iFile},2);
    ax_err(iFile) = 1;
  end
end

err = ax_err;

data = [];
