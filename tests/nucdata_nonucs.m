function ok = test()

Isotopes = '';
[I] = nucdata(Isotopes);
[I,gn] = nucdata(Isotopes);
[I,gn,qm] = nucdata(Isotopes);
[I,gn,qm,ab] = nucdata(Isotopes);

ok = isempty(I) && isempty(gn) && isempty(qm) && isempty(ab);
