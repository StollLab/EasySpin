function ok = test()

BaseDir = 'eprfiles/';
SimpleFile = [BaseDir 'strong1.dta'];

% Check syntax

y = eprload(SimpleFile);
[x,y] = eprload(SimpleFile);
[x,y,p] = eprload(SimpleFile);

ok = true;

