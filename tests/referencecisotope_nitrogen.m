function ok = test()

% Test whether reference isotpes for nitrogen are correct.

[Aref,Qref] = runprivate("referenceisotope","N");

ok = Aref.symbol=="14N" && Qref.symbol=="14N";
