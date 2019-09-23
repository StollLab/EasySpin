function [err,data] = test(opt,olddata)

TestValue = barn;
RealValue = 1e-28;
err = ~areequal(TestValue,RealValue,0,'rel');

data = [];
