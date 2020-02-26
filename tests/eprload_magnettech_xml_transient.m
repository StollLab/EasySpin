function ok = test()

% Read Magnettech spectrometer files (new XML format)
%-------------------------------------------------

BaseDir = 'eprfiles/magnettech/';
Files{1} = 'LOGSTest_006_transient.xml';
fileName = [BaseDir Files{1}];

[x,y,p] = eprload(fileName);

ok = x(1)==0.022 && x(end)==9.976;
