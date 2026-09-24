function ok = test()

% Default number of points must be even, so that the output
% can be passed to rapidscan2spc (issue #260)
%-------------------------------------------------------

g = 2;
T1 = 4;         % µs
T2 = 0.5;       % µs
B1 = 1e-3;      % mT
modAmp = 1;     % mT
modFreq = 50;   % kHz

deltaB0 = [0 0.1];
Methods = {'fft','td'};

ok = true;
for iB = 1:numel(deltaB0)
  for iM = 1:numel(Methods)
    Opt.Method = Methods{iM};
    [t,Mx,My,~] = blochsteady(g,T1,T2,deltaB0(iB),B1,modAmp,modFreq,Opt);
    ok(end+1) = mod(numel(t),2)==0 && numel(Mx)==numel(t) && numel(My)==numel(t);
    M = -Mx + 1i*My;
    [B,spc] = rapidscan2spc(M,modAmp,modFreq,g);
    ok(end+1) = ~isempty(B) && numel(B)==numel(spc);
  end
end
