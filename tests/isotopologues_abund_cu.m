function [err,data] = test(opt,olddata)

ab63 = nucabund('63Cu');
ab65 = nucabund('65Cu');

a = isotopologues('Cu');
b = [ab63 ab65];
ok = areequal(a.Abund,b);

a = isotopologues('Cu,Cu');
b = [ab63^2 ab63*ab65 ab63*ab65 ab65^2];
ok = ok && areequal(a.Abund,b);

a = isotopologues('Cu,Cu,63Cu');
b = [ab63^2 ab63*ab65 ab63*ab65 ab65^2];
ok = ok && areequal(a.Abund,b);

a = isotopologues('Cu,63Cu,Cu');
b = [ab63^2 ab63*ab65 ab63*ab65 ab65^2];
ok = ok && areequal(a.Abund,b);

a = isotopologues('63Cu,Cu,Cu');
b = [ab63^2 ab63*ab65 ab63*ab65 ab65^2];
ok = ok && areequal(a.Abund,b);

a = isotopologues('(63,65)Cu',[.5 .5]);
b = [.5 .5];
ok = ok && areequal(a.Abund,b);

a = isotopologues('(63,65)Cu,1H',{[.5 .5],1});
b = [.5 .5];
ok = ok && areequal(a.Abund,b);

a = isotopologues('(63,65)Cu,(1,2)H',{[.5 .5],[.9 .1]});
b = [.45 .45 .05 .05];
ok = ok && areequal(a.Abund,b);

err = ~ok;
data = [];
