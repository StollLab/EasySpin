function ok = test()

v = runprivate('letter2vec','y');
ok(1) = all(v==[0;1;0]);

str = runprivate('nuclist2string',{'63Cu','1H'});
ok(2) = str=="63Cu,1H";
