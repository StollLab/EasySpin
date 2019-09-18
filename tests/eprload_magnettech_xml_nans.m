function [err,data] = test(opt,olddata)

% Read Magnettech spectrometer files (new XML format)
%-------------------------------------------------

BaseDir = 'eprfiles/magnettech/';

Files{1} = 'LOGSTest_001_field_sweep.xml';

thisFileName = [BaseDir Files{1}];

[x,data] = eprload(thisFileName);

err = any(isnan(x)) || any(isnan(data));

data = [];
