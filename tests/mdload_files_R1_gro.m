function [err,data] = test(opt,olddata)

% Read Gromos89 format
%-------------------------------------------------------------------------------
FileName = 'mdfiles/1ubq_A28R1_noSol';

OutOpt.Verbosity = 0;
OutOpt.keepProtCA = 1;

data = mdload([FileName '.trr'], [FileName '.gro'], [], OutOpt);

ok = data.nSteps==11;
ok = ok & data.dt==2.3e-8;
ok = ok & size(data.dihedrals,1)==5;
err = ~ok;

data = [];

end
