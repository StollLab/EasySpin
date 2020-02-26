function ok = test()

BaseDir = 'eprfiles/';
SimpleFiles = {[BaseDir 'strong1.dta'],[BaseDir 'strong1esp.spc']};

% Check whether the parameter structure is parsed
%-------------------------------------------------
[x,y,p] = eprload(SimpleFiles{1});
ok(1) = isreal(p.XPTS) && isreal(p.XMIN) && isreal(p.XWID);

[x,y,p] = eprload(SimpleFiles{2});
ok(2) = isreal(p.MF) && isreal(p.GST) && isreal(p.MP);
