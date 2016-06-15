function [err,data] = test(opt,olddata)

% Test 2: 
%======================================================
Sz = sop(1,'z');
Sy = sop(1,'y'); % S=1 system
mwFreq = 10e3;
tp = 0.010; % 10 GHz, 10 ns
U = propint(mwFreq*Sz,1/2/tp*Sy,tp,mwFreq);
Ur = [
   0.5000 + 0.0006i  -0.7071 - 0.0004i   0.5000 + 0.0000i
   0.7071 + 0.0004i   0.0000 - 0.0000i  -0.7071 + 0.0004i
   0.5000 + 0.0000i   0.7071 - 0.0004i   0.5000 - 0.0006i
];
err = ~areequal(U,Ur,1e-4);
data = [];
