function ok = test()

% Principal axes x, y and z

thr = 1e-12;

d{1}  = {'x' [1;0;0]};
d{2}  = {'y' [0;1;0]};
d{3}  = {'z' [0;0;1]};
d{4}  = {'xy' [1;1;0]};
d{5}  = {'yz' [0;1;1]};
d{6}  = {'xz' [1;0;1]};
d{7}  = {'xyz' [1;1;1]};

for k = 1:numel(d)
  ang1 = vec2ang(d{k}{1});
  ang2 = vec2ang(d{k}{2});
  ok(k) = areequal(ang1,ang2,thr,'abs');
end
