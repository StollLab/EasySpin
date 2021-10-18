function ok = test()

Isotopes = '112Sn';
[s] = nucdata(Isotopes);
[s,gn] = nucdata(Isotopes);
[s,gn,qm] = nucdata(Isotopes);
[s,gn,qm,ab] = nucdata(Isotopes);

ok = true;
