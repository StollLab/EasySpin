function [err,data] = test(opt,olddata)

% System ------------------------
Sys.ZeemanFreq = 33.500;

% Experiment -------------------
Pulse.tp = 0.03;
Pulse.Flip = pi;

Exp.Sequence = {Pulse 0};
Exp.mwFreq = 33.5;
Exp.DetSequence = [0 1]; 

% Options ---------------------------
Exp.DetOperator = {'z1'};


% one single acquistion point ------------------
[sig] = spidyan(Sys,Exp);

sigexpected = 0.9975441;


% two dimensions -----------------------
Exp.nPoints = [2 3];
Exp.Dim1 = {'p1.tp' 0.005};
Exp.Dim2 = {'p1.Flip' 0.005};

[sig2] = spidyan(Sys,Exp);

data.sig2 = sig2;

% two detection operators ----------------------
Exp.DetOperator = {'+1' 'z1'};
[sig3] = spidyan(Sys,Exp);

data.sig3 = sig3;

if ~isempty(olddata)
  err = [~areequal(sig,sigexpected,1e-4) ~areequal(sig2,olddata.sig2,1e-4)...
    ~areequal(sig3,olddata.sig3,1e-4)];
else
  err = [];
end

