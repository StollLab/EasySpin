function ok = test()

% Test whether reference isotpes for hydrogen are correct.

[Aref,Qref] = runprivate("referenceisotope","H");

ok = Aref.symbol=="1H" && Qref.symbol=="2H";
