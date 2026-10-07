function ok = test()

% Full matrices cannot be combined with nonzero tilt angles (xFrame)

M = [1 0.1 0.05; 0.1 2 0.2; 0.05 0.2 3];  % symmetric full matrix
ang = [0.3 0.7 -0.4];

Sys = {};
Sys{end+1} = struct('S',1/2,'g',2*eye(3)+M/10);
Sys{end+1} = struct('S',1,'D',M*100);
Sys{end+1} = struct('S',[1/2 1/2],'ee',M*10);
Sys{end+1} = struct('S',1/2,'Nucs','1H','A',M*10);
Sys{end+1} = struct('S',1/2,'Nucs','14N','A',[1 1 1],'Q',M-2*eye(3));
Sys{end+1} = struct('S',1/2,'Nucs','1H','A',[1 1 1],'sigma',eye(3)+M*1e-3);
Sys{end+1} = struct('S',1/2,'Nucs','1H,1H','A',[1 1 1;2 2 2],'nn',M);
Fields = {'g','D','ee','A','Q','sigma','nn'};

for k = 1:numel(Sys)
  FrameField = [Fields{k} 'Frame'];

  % full matrix with zero Frame: accepted
  Sys_ = Sys{k};
  Sys_.(FrameField) = [0 0 0];
  [~,err] = runprivate('validatespinsys',Sys_);
  ok(k,1) = isempty(err);

  % full matrix with nonzero Frame: rejected
  Sys_.(FrameField) = ang;
  [~,err] = runprivate('validatespinsys',Sys_);
  ok(k,2) = contains(err,FrameField);
end

% eeFrame is still allowed with J/dip/dvec (internally full ee matrix)
Sys_ = struct('S',[1/2 1/2],'J',10,'dip',[1 2],'dvec',[1 2 3],'eeFrame',ang);
[~,err] = runprivate('validatespinsys',Sys_);
ok(end+1,:) = isempty(err);

ok = ok(:).';
