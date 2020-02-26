function ok = test()

isoList = isotopologues('(63,65)Cu,1H',[2,1],{[0.5 0.5],1});
w = [isoList.weight];
w0 = [0.25 0.5 0.25];
ok = areequal(w,w0,0);
