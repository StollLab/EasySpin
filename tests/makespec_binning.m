function ok = test()

% Peaks go to the nearest point ('round') or the next point below ('floor')

rng(1);
nPoints = 101;
Amp = 4;
Pos = rand(1,100);

ok = true(2,numel(Pos));
for k = 1:numel(Pos)
  [x,s] = makespec([0 1],nPoints,Pos(k),Amp);
  [~,idxNearest] = min(abs(x-Pos(k)));
  ok(1,k) = isequal(find(s),idxNearest) && s(idxNearest)==Amp;

  s = makespec([0 1],nPoints,Pos(k),Amp,'floor');
  idxBelow = find(x<=Pos(k),1,'last');
  ok(2,k) = isequal(find(s),idxBelow) && s(idxBelow)==Amp;
end
ok = all(ok,2);
