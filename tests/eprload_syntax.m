function [err,data] = test(opt,olddata)

BaseDir = 'eprfiles/';
SimpleFile = [BaseDir 'strong1.dta'];

%: Check syntax
%-------------------------------------------------
%eprload(SimpleFile);
[y] = eprload(SimpleFile);
[x y] = eprload(SimpleFile);
[x y p] = eprload(SimpleFile);
err = 0;
data = [];
