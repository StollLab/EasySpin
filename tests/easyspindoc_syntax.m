function ok = test()

docEntry = easyspin('doc');
ref = "index.html";
ok = docEntry(end-strlength(ref)+1:end)==ref;
