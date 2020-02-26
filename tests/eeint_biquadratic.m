function ok = test()

% isotropic biquadratic exchange

j = 10;
maxS = 3;

Sys.ee = 0;
Sys.ee2 = j;

iPair = 1;
for S1 = 0.5:0.5:maxS
  for S2 = 0.5:0.5:maxS
    
    % using EasySpin's sham
    Sys.S = [S1 S2];
    H0 = eeint(Sys);
    
    % manual construction
    [S1x,S1y,S1z] = sop([S1 S2],'xe','ye','ze');
    [S2x,S2y,S2z] = sop([S1 S2],'ex','ey','ez');
    H1 = j*(S1x*S2x + S1y*S2y + S1z*S2z)^2;
    
    ok(iPair) = areequal(H0,H1,1e-10,'abs');
    
    iPair = iPair + 1;
  end
end
