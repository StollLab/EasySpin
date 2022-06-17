function ok = test()

% Test several direction combinations of light polarization and
% molecular orientation

% see
% M.Thurnauer, J. Norris, Biochem. Biophys. Res. Commun. 73, 1976, 501-506
% Figure 1 (only orientations with * give intensity)
% https://doi.org/10.1016/0006-291X(76)90735-X
% (note that they use a left-hand molecular frame!)

% Transition dipole moment
tdm = 'y';  % TDM along molecular y axis

% Light excitation information
k = 'y'; % propagation direction along y
alpha{1} = pi/2; % E-field perpendicular to B0
alpha{2} = 0; % E-field parallel to B0

ori = {[0 0 pi/2], [0 0 0]; [pi/2 pi/2 0], [0 pi/2 0]; [pi/2 pi/2 pi/2], [0 pi/2 pi/2]}; 

for iDir = 1:2
  for iOri1 = 1:3
    for iOri2 = 1:2
      ori_ = -ori{iOri1,iOri2}(3:-1:1);
      weights{iDir}(iOri1,iOri2) = photoselect(tdm,ori_,k,alpha{iDir});
    end
  end
end

weightsRef{1} = [1 0; 1 0; 0 0];
weightsRef{2} = [0 0; 0 0; 1 1];

threshold = 1e-10;
ok(1) = areequal(weights{1},weightsRef{1},threshold,'abs');
ok(2) = areequal(weights{2},weightsRef{2},threshold,'abs');
