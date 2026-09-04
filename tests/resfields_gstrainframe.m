function ok = test

Sys.S = 1/2;
Sys.g = [2.000 2.005 2.010];
Sys.gStrain = [0.0001 0.0003 0.005];

Exp.mwFreq = 34; % GHz
Exp.Range = [1204 1217]; % mT
Exp.Harmonic = 0;

gFrame = [0 0 0; -20 90 0; -20 50 80; -90 30 0]*pi/180;
nFrames = size(gFrame,1);

for i = 1:nFrames
  Sys.gFrame = gFrame(i,:);
  Exp.SampleFrame = -fliplr(gFrame(i,:));
  [B(i),I(i),W(i)] = resfields(Sys,Exp);
end

for k = 1:nFrames-1
  ok(k) = areequal(W(k),W(end),1e-8,'rel');
end

end
