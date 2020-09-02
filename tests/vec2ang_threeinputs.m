function ok = test()

% Test three inputs x,y,z versus one input v.

n = 1000;
x = rand(1,n);
y = rand(1,n);
z = rand(1,n);

v = [x;y;z];

a = vec2ang(v);
b = vec2ang(x,y,z);

ok = areequal(a,b,1e-10,'abs');
