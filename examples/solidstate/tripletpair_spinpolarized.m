% Polarization of a triplet pair generated through singlet fission 
%=================================================================
%
% Parameters adapted from
% Tayebjee, M., Sanders, S., Kumarasamy, E. et al. 
% Nature Phys 13, 182–188 (2017). 
% https://doi.org/10.1038/nphys3909

clear

Sys.S = [1 1];
Sys.g = [2.002 2.002];
Sys.lwpp =1;

Sys.D = [1140 20;1140 20]; 
Sys.J = 20000;

% [++ +0 +- 0+ 00 0- -+ -0 --]
Sys.Pop = [0 0 1 0 1 0 1 0 0];
% With Sys.PopBasis = 'Spin', the populations of the sublevels are 
% calculated as the overlap between the eigenvectors and the provided state vector  
Sys.PopMode = 'eigenbasis'; 

Exp.Harmonic = 0;
Exp.mwFreq = 9.5;%
Exp.Range = [290 390];
Exp.nPoints = 1024;

Opt = [];

figure(1)
pepper(Sys,Exp,Opt);
figure(2)
ori = 'xyz';

B = 290:390;
E = zeros(length(B),9);
for i = 1:3
  subplot(3,1,i)
  levelsplot(Sys,ori(i),B,9.5)
  title(ori(i))
end
hold off

