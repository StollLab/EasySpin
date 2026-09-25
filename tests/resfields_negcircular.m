function ok = test()

% Circular polarization at negative fields: the sense of rotation relative
% to the field is reversed

Sys.g = 2;
Exp.mwFreq = 9.5;
Exp.SampleFrame = [0.3 0.4 0.5];

Exp.Range = [330 350];
Exp.mwMode = {0 'circular+'};
[~,Ipp] = resfields(Sys,Exp);
Exp.mwMode = {0 'circular-'};
[~,Ipm] = resfields(Sys,Exp);

Exp.Range = [-350 -330];
Exp.mwMode = {0 'circular+'};
[~,Inp] = resfields(Sys,Exp);
Exp.mwMode = {0 'circular-'};
[~,Inm] = resfields(Sys,Exp);

Int = [Ipp Ipm Inp Inm]/Ipp;

ok = areequal(Int,[1 0 0 1],1e-6,'abs');
