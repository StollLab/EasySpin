function ok = test()

% Test 1: syntax
%======================================================
x = linspace(0,10,201);
y = 3*exp(-x);
[k] = exponfit(x,y);
[k,c] = exponfit(x,y);
[k,c,yFit] = exponfit(x,y);
[k,c] = exponfit(x,y,1);
[k,c] = exponfit(x,y,1,'noconst');

ok = true;

