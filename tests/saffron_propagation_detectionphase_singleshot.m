function [ok,data] = test(opt,olddata) %#ok<FNDEF>

% Detection phase: pi negates the signal of the default detection operator,
% and with several operators each one has its own phase (single
% orientation, detection over the whole sequence)

Sys.S = 1/2;

p.Flip = pi/2;
p.tp = 0.02;

Exp.Sequence = {p 0.5};
Exp.mwFreq = 9.5;
Exp.Field = mhz2mt(9500);  % mT, resonance field for g = gfree
Exp.MolFrame = [0 0 0];  % single orientation
Exp.DetWindow = 'sequence';

% Test 1 - comparison for default detection operator
s1_nodetectionphase = saffron(Sys,Exp);

Exp_ = Exp;
Exp_.DetPhase = pi;

s1_detectionphase = saffron(Sys,Exp_);

if (opt.Display)
    subplot(3,1,[1 2]);
    hold on
    plot(real(-s1_nodetectionphase),'r')
    plot(real(s1_detectionphase),'b')
    axis tight
    legend('no detection phase','with detection phase');
    title([mfilename '- default detection operator']);
    subplot(3,1,3);
    plot(real(s1_nodetectionphase+s1_detectionphase));
    axis tight
    xlabel('data points direct dimension');
end

% Test 2 - several detection operators

Exp_ = Exp;

Exp_.DetOperator = {'z1' '+1' '+1'};
Exp_.DetPhase = [0 0 pi];

s2 = saffron(Sys,Exp_);

data.s2 = s2;

% Comparison

if ~isempty(olddata)
  ok = [areequal(s1_detectionphase,-s1_nodetectionphase,1e-4,'abs') ... % test 1
        areequal(s2(:),olddata.s2(:),1e-4,'abs') areequal(s2(:,2),-s2(:,3),1e-4,'abs')]; % test 2
else
  ok = [];
end
