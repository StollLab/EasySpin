function ok = test(opt)

% Check for correct conversion of [px py pz] to energy-ordered zero-field
% populations for given signs of D and E

DE(1,:) = [ 900 -200];  % MHz
DE(2,:) = [-600   50];  % MHz
DE(3,:) = [ 800  100];  % MHz
DE(4,:) = [-200  -55];  % MHz

for i = 1:size(DE,1)
  D = DE(i,1);
  E = DE(i,2);

  Triplet.S = 1;
  Triplet.D = [D E];  % MHz

  % Provide population vector as [px py pz]
  xyzpops = rand(1,3);

  % Determine energy-ordered populations using D and E
  popsDE = runprivate('tripletpoporder',D,E,xyzpops);

  % Determine energy-ordered populations using eigenvalues for ZFS
  % interaction (conversion happens in validatespinsys)
  Sys = runprivate('validatespinsys',Triplet);
  popsD = runprivate('tripletpoporder',Sys.D,xyzpops);

  % Corresponding energy-ordered population vector
  px = xyzpops(1); py = xyzpops(2); pz = xyzpops(3);
  if D>0 && E<0
    p1 = pz; p2 = py; p3 = px;
  elseif D>0 && E>0
    p1 = pz; p2 = px; p3 = py;
  elseif D<0 && E>0
    p1 = px; p2 = py; p3 = pz;
  elseif D<0 && E<0
    p1 = py; p2 = px; p3 = pz;
  end
  pops0 = [p1 p2 p3];

  ok(i) = areequal(pops0,popsDE,1e-16,'abs') && areequal(pops0,popsD,1e-16,'abs');

end

end
