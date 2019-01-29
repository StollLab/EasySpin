function [err,data] = test(opt,olddata)

Sys.ZeemanFreq = 9.5;

p.Flip = pi/2;
p.tp = 0.02;

Exp.Sequence = {p 0.5};
Exp.mwFreq = 9.5;

%% Test 1 - comparison for default detection operator
s1_nodetectionphase = spidyan(Sys,Exp);

Exp_ = Exp;
Exp_.DetPhase = pi;

s1_detectionphase = spidyan(Sys,Exp_);

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

%% Test 2 - several detection operators

Exp_ = Exp;

Exp_.DetOperator = {'z1' '+1' '+1'};
Exp_.DetFreq = [0 9.5 9.5];
Exp_.DetPhase = [0 0 pi];

s2 = spidyan(Sys,Exp_);

data.s2 = s2;

%% Comparison

if ~isempty(olddata)
  err = [~areequal(s1_detectionphase,-s1_nodetectionphase,1e-4) ... % comparison of Test1
         ~areequal(s2,olddata.s2,1e-4) ~areequal(s2(:,2),-s2(:,3),1e-4)]; % test 2
else
  err = [];
end

