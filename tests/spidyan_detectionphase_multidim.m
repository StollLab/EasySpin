function ok = test()

Sys.ZeemanFreq = 9.5;

p.Flip = pi/2;
p.tp = 0.02;

Exp.Sequence = {p 0.5};
Exp.mwFreq = 9.5;
Exp.DetOperator = {'z1' '+1'};
Exp.DetFreq = [0 9.5];

reference = spidyan(Sys,Exp);

% Test if phase is correctly applied to numeric arrays
Exp_ = Exp;
Exp_.nPoints = 2;
Exp_.Dim1 = {'p1.Flip' pi/4};
Exp_.DetPhase = [pi pi];

test1 = spidyan(Sys,Exp_);

% Test if phase is correctly applied to cell arrays
Exp_ = Exp;
Exp_.nPoints = 2;
Exp_.Sequence{2} = 0.3;
Exp_.Dim1 = {'d1' 0.2};
Exp_.DetPhase = [pi pi];

test2 = spidyan(Sys,Exp_);


% Comparison
ok = areequal(reference,-squeeze(test1(1,:,:)),1e-4,'abs') &&...
     areequal(reference,-test2{2},1e-4,'abs');
