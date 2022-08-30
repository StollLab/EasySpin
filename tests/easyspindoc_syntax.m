function ok = test()

docEntry = easyspindoc;
ref = "index.html";
ok = docEntry(end-strlength(ref)+1:end)==ref;
