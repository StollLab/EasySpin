function [err,data] = test(opt,olddata)

Freq = 1500;
DetectionOperators = {'+1'};

t = 0:0.0001:0.2;
TimeAxis{1} = t;
Cos = cos(2*pi*Freq*t);
Sin = sin(2*pi*Freq*t);

RawSignal1{1} = Cos;

RawSignal2{1} = [Cos; Cos + 1i*Sin; ones(1,length(t)); ones(1,length(t)); ones(1,length(t))];
RawSignal2{2} = [Cos + 1i*Sin*1e-10; Cos*1e-10 + 1i*Sin; ones(1,length(t)); ones(1,length(t))];


[ProcessedSignal1] = signalprocessing(TimeAxis,RawSignal1,DetectionOperators,-Freq/1000);

data.ProcessedSignal1 = ProcessedSignal1;
TimeAxis{2} = t;

DetectionOperators = {'+1','x2','x2','x2'};
[ProcessedSignal2] = signalprocessing(TimeAxis,RawSignal2,DetectionOperators,[-Freq/1000 -Freq/1000 0]);

data.ProcessedSignal2 = ProcessedSignal2;


if ~isempty(olddata)
  err = ~areequal(ProcessedSignal1{1},olddata.ProcessedSignal1{1},1e-4);
  for i = 1 : length(ProcessedSignal2)
    err(i+1) = ~areequal(ProcessedSignal2{i},olddata.ProcessedSignal2{i},1e-4);
  end
else
  err = [];
end

