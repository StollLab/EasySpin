function [err,data] = test(opt,olddata)

% Compare the two slightly different version of the SPE
% binary format from Magnettech spectrometers
%---------------------------------------------------------

[d1,d2,par1] = eprload('eprfiles/magnettech/oldformat.spe');
[d1,d2,par2] = eprload('eprfiles/magnettech/newformat.spe');

ok = areequal(par1.B0_Field,par2.B0_Field,1e-3);
ok = ok && areequal(par1.B0_Scan,par2.B0_Scan,1e-3);
err = ~ok;

data = [];
