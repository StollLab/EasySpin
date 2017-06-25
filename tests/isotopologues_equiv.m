function [err,data] = test(opt,olddata)

isoList = isotopologues('(63,65)Cu,1H',[2,1],{[0.5 0.5],1});
abund = extractabundances(isoList);
abund0 = [0.25 0.5 0.25];
ok = areequal(abund,abund0);

err = ~ok;
data = [];

function a = extractabundances(iso)
a = [iso.Abund];
