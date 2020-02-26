function ok = test()

% simple fits with and without offset
%======================================================
x = linspace(0,1,301);
y = 7*exp(-10*x) + 2*exp(-1*x) + 0.2*rand(size(x));

[k,c] = exponfit(x,y,1);
[K,C] = exponfit(x,y,1,'noconst');
[k,c] = exponfit(x,y,2);
[K,C] = exponfit(x,y,2,'noconst');

ok = true;
