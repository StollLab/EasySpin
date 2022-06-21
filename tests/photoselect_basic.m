function ok = test()

% Test photoselection values for all direction combinations of light polarization
% and molecular orientation against values from literature.

% See
% M.Thurnauer, J. Norris, Biochem. Biophys. Res. Commun. 73, 1976, 501-506
% https://doi.org/10.1016/0006-291X(76)90735-X
% Figure 1 (only orientations with * give intensity)
% (note that they use a left-hand molecular frame!)

% Transition dipole moment
tdm = 'y';  % tdm along molecular y axis

% Light excitation information
k = 'y'; % propagation direction along laboratory y axis
alpha{1} = pi/2; % E-field along laboratory x axis (perpendicular to B0)
alpha{2} = 0; % E-field along laboratory z axis (parallel to B0)

% List of molecular orientation as in the figure
ori = {[0 0 pi/2], [0 0 0]; [pi/2 pi/2 0], [0 pi/2 0]; [pi/2 pi/2 pi/2], [0 pi/2 pi/2]}; 
ori = cellfun(@(x)-x(3:-1:1),ori,'UniformOutput',false);

for iDir = 1:2
  for iOri = 1:3
    for jOri = 1:2
      weights{iDir}(iOri,jOri) = photoselect(tdm,ori{iOri,jOri},k,alpha{iDir});
    end
  end
end

weights_reference{1} = [1 0; 1 0; 0 0];
weights_reference{2} = [0 0; 0 0; 1 1];

threshold = 1e-10;
ok(1) = areequal(weights{1},weights_reference{1},threshold,'abs');
ok(2) = areequal(weights{2},weights_reference{2},threshold,'abs');
