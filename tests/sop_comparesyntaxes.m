function [err,data] = test(opt,olddata)

% Compare the two different character-based syntaxes 
%
%  Syntax 1:  'eex'
%  Syntax 2:  'x3'

testfun = @(spins,syntax1,syntax2) areequal(sop(spins,syntax1),sop(spins,syntax2),1e-10);

ok = true;
ok = ok && testfun(1,'x','x1');
ok = ok && testfun([1/2 1],'xe','x1');
ok = ok && testfun([1/2 1],'xy','x1y2');
ok = ok && testfun([1/2 1],'e+','+2');
ok = ok && testfun([1/2 1 1/2],'ezz','z2z3');
ok = ok && testfun([1/2 1 1/2],'ezz','z3z2');

err = ~ok;

data = [];
