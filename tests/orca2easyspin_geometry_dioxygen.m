function ok = test()

% Check that tensor principal axes are oriented correctly relative to the
% molecular structure, for all supported file types:
%   dioxygen: unique axes of g and 17O A tensors along O-O bond

files = {
  'v4.2.1/dioxygen.out'
  'v4.2.1/dioxygen.prop'
  'v5.0.4/dioxygen.out'
  'v6.1.0/dioxygen.out'
  'v6.1.0/dioxygen.property.txt'
  };

% i-th principal axis of tensor with frame ang, in molecular frame
paxis = @(ang,i) erot(ang).'*((1:3).'==i);
% index of unique principal value of an axial tensor
iunique = @(v) find(abs(v-median(v))==max(abs(v-median(v))),1);
% absolute cosine between two vectors
cosabs = @(u,v) abs(u.'*v)/norm(u)/norm(v);

for k = numel(files):-1:1
  Sys = orca2easyspin(['orca/' files{k}]);
  bond = (Sys.xyz(2,:)-Sys.xyz(1,:)).';
  gz = paxis(Sys.gFrame,iunique(Sys.g));
  okk = cosabs(gz,bond)>1-1e-6;
  for n = 1:2
    Az = paxis(Sys.AFrame(n,:),iunique(Sys.A(n,:)));
    okk = okk && cosabs(Az,bond)>1-1e-6;
  end
  ok(k) = okk;
end
