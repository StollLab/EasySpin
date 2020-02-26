function ok = test()

% 3x S=1/2, full ee tensors given
% (Just should not crash, as it did in version 2.6.0)

Sys = struct('S',[1/2 1/2 1/2],'g',[2 2 2]);

J1=-100; %1-2, 2-3
J2=-89.9; % 1-3
dz1=4.85;
dz2=-dz1;

ee12 = [-J1 dz1 0; -dz1 -J1 0; 0 0 -J1];
ee23 = ee12;
ee13 = [-J2 dz2 0; -dz2 -J2 0; 0 0 -J2];
EE = [ee12;ee13;ee23]; % three coupling matrices
Sys.ee=EE*100*clight/1e6; % conversion cm^-1 -> MHz

G= symm(Sys);

ok = true;
