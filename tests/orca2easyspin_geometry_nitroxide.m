function ok = test()

% Check that tensor principal axes are oriented correctly relative to the
% molecular structure, for all supported file types:
%   nitroxide: largest 14N A principal axis and smallest g principal axis
%     perpendicular to C-N-C plane, largest g principal axis along N-O bond

files = {
  'v4.2.1/nitroxide.out'
  'v4.2.1/nitroxide.prop'
  'v5.0.4/nitroxide.out'
  'v6.1.0/nitroxide.out'
  'v6.1.0/nitroxide.property.txt'
  };

% i-th principal axis of tensor with frame ang, in molecular frame
paxis = @(ang,i) erot(ang).'*((1:3).'==i);
% absolute cosine between two vectors
cosabs = @(u,v) abs(u.'*v)/norm(u)/norm(v);

for k = numel(files):-1:1
  Sys = orca2easyspin(['orca/' files{k}]);
  C3 = Sys.xyz(3,:); N = Sys.xyz(4,:); C5 = Sys.xyz(5,:); O = Sys.xyz(6,:);
  normal = cross(C3-N,C5-N).';
  NO = (O-N).';
  iN = find(Sys.NucsIdx==4);
  [~,iAz] = max(abs(Sys.A(iN,:)));
  [~,igz] = min(Sys.g);
  [~,igx] = max(Sys.g);
  okk = cosabs(paxis(Sys.AFrame(iN,:),iAz),normal)>0.98;
  okk = okk && cosabs(paxis(Sys.gFrame,igz),normal)>0.98;
  okk = okk && cosabs(paxis(Sys.gFrame,igx),NO)>0.98;
  ok(k) = okk;
end
