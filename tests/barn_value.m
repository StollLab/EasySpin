function [err,data] = test(opt,olddata)

TestValue = barn;
RealValue = 1e-28;
err = ~areequal(TestValue,RealValue);

data = [];
