function [err,data] = test(opt,olddata)

Freq = 1500;

TimeAxis = 0:0.0001:0.2;
Cos = cos(2*pi*Freq*TimeAxis);
Sin = sin(2*pi*Freq*TimeAxis);

RawSignal1 = Cos.';

RawSignal2(1,:,:) = [Cos; Cos + 1i*Sin; ones(1,length(TimeAxis)); ones(1,length(TimeAxis))];
RawSignal2(2,1:2,:) = [Cos + 1i*Sin*1e-10; Cos*1e-10 + 1i*Sin];

[ProcessedSignal1] = signalprocessing(TimeAxis,RawSignal1,-Freq/1000);
[ProcessedSignal1_] = signalprocessing(TimeAxis,RawSignal1.',-Freq/1000);

data.ProcessedSignal1 = ProcessedSignal1;

[ProcessedSignal2] = signalprocessing(TimeAxis,RawSignal2,[-Freq/1000 -Freq/1000 0 0]);

data.ProcessedSignal2 = ProcessedSignal2;


if ~isempty(olddata)
  err = [~isequal(ProcessedSignal1,ProcessedSignal1_) ~areequal(ProcessedSignal1,olddata.ProcessedSignal1,1e-4,'abs') ~areequal(ProcessedSignal2,olddata.ProcessedSignal2,1e-4,'abs')];
else
  err = [];
end

