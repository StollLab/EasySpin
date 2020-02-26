function ok = test()

ab63 = nucabund('63Cu');
ab65 = nucabund('65Cu');

extractweights = @(iso)[iso.weight];

thr = 1e-10;

isoList = isotopologues('Cu');
w0 = [ab63 ab65];
w = extractweights(isoList);
ok = areequal(w,w0,thr,'abs');

isoList = isotopologues('Cu,Cu');
w = extractweights(isoList);
w0 = [ab63^2 ab63*ab65 ab63*ab65 ab65^2];
ok = ok && areequal(w,w0,thr,'abs');

isoList = isotopologues('Cu,Cu,63Cu');
w = extractweights(isoList);
w0 = [ab63^2 ab63*ab65 ab63*ab65 ab65^2];
ok = ok && areequal(w,w0,thr,'abs');

isoList = isotopologues('Cu,63Cu,Cu');
w = extractweights(isoList);
w0 = [ab63^2 ab63*ab65 ab63*ab65 ab65^2];
ok = ok && areequal(w,w0,thr,'abs');

isoList = isotopologues('63Cu,Cu,Cu');
w = extractweights(isoList);
w0 = [ab63^2 ab63*ab65 ab63*ab65 ab65^2];
ok = ok && areequal(w,w0,thr,'abs');

isoList = isotopologues('(63,65)Cu',1,[0.5 0.5]);
w = extractweights(isoList);
w0 = [0.5 0.5];
ok = ok && areequal(w,w0,thr,'abs');

isoList = isotopologues('(63,65)Cu,1H',[1 1],{[0.5 0.5],1});
w = extractweights(isoList);
w0 = [0.5 0.5];
ok = ok && areequal(w,w0,thr,'abs');

isoList = isotopologues('(63,65)Cu,(1,2)H',[1,1],{[0.5 0.5],[0.9 0.1]});
w = extractweights(isoList);
w0 = [0.45 0.05 0.45 0.05];
ok = ok && areequal(w,w0,thr,'abs');
