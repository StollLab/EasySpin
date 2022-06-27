function ok = test()

% photoselect() should be able to handle multiple orientations

rng(1);

tdm = [23 34]*pi/180;
k = [155 244]*pi/180;
alpha = 12*pi/180;

nOrientations = 7;
ori = rand(nOrientations,3).*[2 1 2]*pi;

w = photoselect(tdm,ori,k,alpha);

ok = numel(w)==nOrientations;
