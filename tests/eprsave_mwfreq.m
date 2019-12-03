function [err,data] = test(opt,olddata)

%======================================================
% Save and load 1D dataset with microwave frequency
%======================================================

fn = tempname;
mwFreq = 9.12345;

N = 1001;
x = linspace(300,400,N);
y = rand(1,N);
eprsave(fn,x,y,'title',mwFreq);
[xo,yo,p] = eprload(fn);
delete([fn '.*']);

mwFreqLoaded = p.MWFQ/1e9;

err = ~areequal(mwFreq,mwFreqLoaded,1e-9,'abs');

data = [];
