function [ok,data] = test(opt,olddata) %#ok<FNDEF>

% Single-point detection at the end of the sequence, without and with two
% indirect dimensions, and with two detection operators (single
% orientation)

% System ------------------------
Sys.S = 1/2;

% Experiment -------------------
Pulse.tp = 0.03;
Pulse.Flip = pi;

Exp.Sequence = {Pulse};
Exp.mwFreq = 33.5;
Exp.Field = mhz2mt(33500);  % mT, resonance field for g = gfree
Exp.MolFrame = [0 0 0];  % single orientation
Exp.DetWindow = 0;

% Options ---------------------------
Exp.DetOperator = {'z1'};


% one single acquistion point ------------------
sig = saffron(Sys,Exp);

sigexpected = 1.0000;


% two dimensions -----------------------
Exp.nPoints = [2 3];
Exp.Dim1 = {'p1.tp' 0.005};
Exp.Dim2 = {'p1.Flip' 0.005};

sig2 = saffron(Sys,Exp);

data.sig2 = sig2;

% two detection operators ----------------------
Exp.DetOperator = {'+1' 'z1'};
sig3 = saffron(Sys,Exp);

data.sig3 = sig3;

if ~isempty(olddata)
  ok = [areequal(sig,sigexpected,1e-4,'abs') areequal(sig2(:),olddata.sig2(:),1e-4,'abs') ...
    areequal(sig3(:),olddata.sig3(:),1e-4,'abs')];
else
  ok = [];
end
