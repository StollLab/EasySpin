function ok = test() %#ok<FNDEF>

% Detection phases are applied to signals with an indirect dimension: a
% detection phase of pi negates S+ and leaves Sz unchanged (single
% orientation, detection over the whole sequence)

Sys.S = 1/2;

p.Flip = pi/2;
p.tp = 0.02;

Exp.Sequence = {p 0.5};
Exp.mwFreq = 9.5;
Exp.Field = mhz2mt(9500);  % mT, resonance field for g = gfree
Exp.MolFrame = [0 0 0];  % single orientation
Exp.DetWindow = 'sequence';
Exp.DetOperator = {'z1' '+1'};

reference = saffron(Sys,Exp);

Exp_ = Exp;
Exp_.nPoints = 2;
Exp_.Dim1 = {'p1.Flip' pi/4};
Exp_.DetPhase = [pi pi];

test1 = saffron(Sys,Exp_);

flip = [1 -1];
ok = areequal(reference.*flip,squeeze(test1(1,:,:)),1e-4,'abs');
