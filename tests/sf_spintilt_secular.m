function ok = test() %#ok<FNDEF>

% Tilt of the electron spin by a weak transverse field (sf_spintilt), for
% the hyperfine enhancement of RF pulses.
% - S = 1/2, anisotropic g: <a|dS|a> is the derivative of the spin
%   expectation value of each level with respect to the field b (finite
%   differences), and equals mS*(g.'*n - q*(q.'*g.'*n))/(|g.'*zL|*B).
% - Secular states (dSsec), S = 1/2 with anisotropic g and S = 5/2 with
%   zero-field splitting: no elements between different mS, and the
%   diagonal holds <a|dS|a> of the eigenstate assigned to each mS. dS
%   itself has a smaller diagonal (cos(theta)*<a|dS|a> for S = 1/2).

B = 350;  % mT
n = [1;0;0];  % direction of the weak field (x of the frame of Sys)

% S = 1/2, anisotropic g
Sys.S = 1/2;
Sys.g = [1.9 2.1 2.6];
Sys.gFrame = [0 0.8 0];
[H0,mux,~,muz] = ham(Sys);
H = full(H0 - B*muz);
V = -full(mux);
Sops = {full(sop(Sys.S,'x')),full(sop(Sys.S,'y')),full(sop(Sys.S,'z'))};
[dS,dSsec] = runprivate('sf_spintilt',H,V,Sops,Sops{3});
[U,~] = eig(H,'vector');

% finite differences of <S> of each level (levels ordered by energy)
db = 1e-4;  % mT
[Up,~] = eig(H+db*V,'vector');
[Um,~] = eig(H-db*V,'vector');
dSnum = zeros(3,2); dSa = zeros(3,2);
for a = 1:2
  for k = 1:3
    dSnum(k,a) = real(Up(:,a)'*Sops{k}*Up(:,a) - Um(:,a)'*Sops{k}*Um(:,a))/(2*db);
    dSa(k,a) = real(U(:,a)'*dS{k}*U(:,a));
  end
end
ok(1) = areequal(dSa,dSnum,1e-6*max(abs(dSnum(:))),'abs');

% analytic tilt of the quantization axis
R = erot(Sys.gFrame).';
g = R*diag(Sys.g)*R.';
zL = [0;0;1];
q = g.'*zL/norm(g.'*zL);
dq = (g.'*n - q*(q.'*g.'*n))/(norm(g.'*zL)*B);
mS = [-1/2 1/2];  % lower level: -q/2 (positive g)
ok(2) = areequal(dSa,dq*mS,1e-9*max(abs(dq)),'abs');

% secular states: diagonal in Sz, with the eigenstate values
ok(3) = checkSecular(dS,dSsec,U,Sops{3});
ok(4) = abs(real(dS{1}(1,1)))<0.999*abs(real(dSsec{1}(1,1)));  % dS is diluted

% S = 5/2 with zero-field splitting
SysM.S = 5/2;
SysM.g = 2;
SysM.D = [600 100];  % MHz
SysM.DFrame = [0 0.3 0];
[H0,mux,~,muz] = ham(SysM);
H = full(H0 - B*muz);
SopsM = {full(sop(SysM.S,'x')),full(sop(SysM.S,'y')),full(sop(SysM.S,'z'))};
[dS,dSsec] = runprivate('sf_spintilt',H,-full(mux),SopsM,SopsM{3});
[U,~] = eig(H,'vector');
ok(5) = checkSecular(dS,dSsec,U,SopsM{3});

end

% dSsec has no elements between different eigenvalues of Sz (diagonal),
% and its diagonal holds <a|dS|a> of the eigenstate a with the largest
% weight on each Sz basis state
function ok = checkSecular(dS,dSsec,U,Sz)
mS = diag(Sz);
ok = true;
[~,assigned] = max(abs(U).^2,[],2);  % eigenstate assigned to each basis state
scale = max(cellfun(@(M) max(abs(M(:))),dSsec));  % some components can be zero
for k = 1:numel(dS)
  offDiagonal = dSsec{k}.*(abs(mS-mS.')>1e-10);
  ok = ok && max(abs(offDiagonal(:)))<=1e-12*scale;
  for i = 1:numel(mS)
    a = assigned(i);
    ok = ok && abs(dSsec{k}(i,i) - U(:,a)'*dS{k}*U(:,a))<=1e-10*scale;
  end
end
end
