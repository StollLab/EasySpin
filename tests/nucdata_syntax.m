function ok = test()

Isotopes = '112Sn';
[I] = nucdata(Isotopes);
[I,gn] = nucdata(Isotopes);
[I,gn,qm] = nucdata(Isotopes);
[I,gn,qm,ab] = nucdata(Isotopes);

ok = true;
