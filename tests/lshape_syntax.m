function ok = test()

x = linspace(-100,100,1e3);
x0 = 34;
w = 20;

y0 = lshape(x,x0,w);
y1 = lshape(x,x0,w,1);
y2 = lshape(x,x0,w,2);
y3 = lshape(x,x0,w,1,0);
y4 = lshape(x,x0,w,1,0.5);
y5 = lshape(x,x0,w,1,1);

ok = true;
