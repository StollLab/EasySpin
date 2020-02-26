function ok = test(opt)

% Simple test of post-convolution vs full treatment

Sys.g = [2.08 2.006 2.002];
Sys.Nucs = '14N,1H';
Sys.A = [20 20 100; 5 5 8];
Sys.tcorr = 0.02e-9;

Exp.mwFreq = 9.5;
Exp.Range = [330 338];
Exp.nPoints = 1e4;

Opt.LLMK = [6 0 2 2];

Opt.PostConvNucs = [];
[x,y1] = chili(Sys,Exp,Opt);
Opt.PostConvNucs = 2;
[x,y2] = chili(Sys,Exp,Opt);

N = 200;
idx = N:Exp.nPoints-N;  % chop off convolution artifacts
y1 = y1(idx+1); % kludge - unclear why the two spectra are one point offset from each other
y2 = y2(idx);
x = x(idx);

if opt.Display
  plot(x,y1,x,y2)
  legend('SLE only','SLE + postconv');
  axis tight
end

% Calculate maximum relative difference
mrd = max(abs(y2(:)-y1(:)))/max(abs(y1(:)));

ok = mrd<=0.02;
