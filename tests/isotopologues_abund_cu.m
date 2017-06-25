function [err,data] = test(opt,olddata)

ab63 = nucabund('63Cu');
ab65 = nucabund('65Cu');

isoList = isotopologues('Cu');
abund0 = [ab63 ab65];
abund = extractabundances(isoList);
ok = areequal(abund,abund0);

isoList = isotopologues('Cu,Cu');
abund = extractabundances(isoList);
abund0 = [ab63^2 ab63*ab65 ab63*ab65 ab65^2];
ok = ok && areequal(abund,abund0);

isoList = isotopologues('Cu,Cu,63Cu');
abund = extractabundances(isoList);
abund0 = [ab63^2 ab63*ab65 ab63*ab65 ab65^2];
ok = ok && areequal(abund,abund0);

isoList = isotopologues('Cu,63Cu,Cu');
abund = extractabundances(isoList);
abund0 = [ab63^2 ab63*ab65 ab63*ab65 ab65^2];
ok = ok && areequal(abund,abund0);

isoList = isotopologues('63Cu,Cu,Cu');
abund = extractabundances(isoList);
abund0 = [ab63^2 ab63*ab65 ab63*ab65 ab65^2];
ok = ok && areequal(abund,abund0);

isoList = isotopologues('(63,65)Cu',1,[0.5 0.5]);
abund = extractabundances(isoList);
abund0 = [0.5 0.5];
ok = ok && areequal(abund,abund0);

isoList = isotopologues('(63,65)Cu,1H',[1 1],{[0.5 0.5],1});
abund = extractabundances(isoList);
abund0 = [0.5 0.5];
ok = ok && areequal(abund,abund0);

isoList = isotopologues('(63,65)Cu,(1,2)H',[1,1],{[0.5 0.5],[0.9 0.1]});
abund = extractabundances(isoList);
abund0 = [0.45 0.05 0.45 0.05];
ok = ok && areequal(abund,abund0);

err = ~ok;
data = [];

function a = extractabundances(iso)
a = [iso.Abund];
