function [err,data] = test(opt,olddata)

% Compare the two slightly different version of the SPE
% binary format from Magnettech spectrometers
%---------------------------------------------------------

[~,~,par1] = eprload('eprfiles/oldformat.spe');
[~,~,par2] = eprload('eprfiles/newformat.spe');

ok = areequal(par1.B0_Field,par2.B0_Field,1e-3);
ok = ok && areequal(par1.B0_Scan,par2.B0_Scan,1e-3);
err = ~ok;

data = [];
