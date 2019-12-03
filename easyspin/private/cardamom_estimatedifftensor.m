function [Diff, msadp, tLag] = estimatedifftensor(RTraj_L, dt, stopFitT)

nSteps = length(RTraj_L);
          
RAlign = RTraj_L(:,:,1);

RTraj_M = zeros(3,3,nSteps);

for iStep = 1:nSteps
  RTraj_M(:,:,iStep) = RTraj_L(:,:,iStep)*RAlign.';
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

RRot = zeros(3,3,nSteps-1);
for iStep = 1:nSteps-1
  RRot(:,:,iStep) = RTraj_M(:,:,iStep+1)*RTraj_M(:,:,iStep).';
end

qRot = rotmat2quat(RRot);
qRot = cat(2, [1;0;0;0], qRot);

idx = qRot(1,:) < 0;
qRot(:,idx) = -qRot(:,idx);

% calculate Cartesian angular velocity components in the molecular frame
wp = q2wp(qRot, dt);

Deltawp = integral(wp, dt);

msadp = msd_fft(Deltawp);
msadp = msadp(:, 1:round(end/2));

tLag = dt*(0:length(msadp)-1)/1e-9;

stopFitN = ceil(stopFitT/dt);

pxp = polyfit(tLag(1:stopFitN), msadp(1,1:stopFitN), 1);
pyp = polyfit(tLag(1:stopFitN), msadp(2,1:stopFitN), 1);
pzp = polyfit(tLag(1:stopFitN), msadp(3,1:stopFitN), 1);

Diff = [pxp(1), pyp(1), pzp(1)]/2*1e9;

end

% Helper functions
% -------------------------------------------------------------------------
function dy = derivative(y, dt)
  dy = zeros(size(y));
  dy(:,2:end-1) = (y(:,3:end) - y(:,1:end-2));
  dy(:,1) = 4*y(:,2) - 3*y(:,1) - y(:,3);
  dy(:,end) = 3*y(:,end) + y(:,end-2) - 4*y(:,end-1);
  dy = dy./(2*dt);
end

function iy = integral(y, dt)
  iy = zeros(size(y));
  iy(:,1) = 0;
  iy(:,2:end-1) = 5*y(:,1:end-2) + 8*y(:,2:end-1) - y(:,3:end);
  iy(:,end) = -y(:,end-2) + 8*y(:,end-1) + 5*y(:,end);
  iy = cumsum(iy, 2)*dt/12;
end

function w = q2w(qTraj, dt)

dq = derivative(qTraj, dt);

q0 = qTraj(1,:,:);
q1 = qTraj(2,:,:);
q2 = qTraj(3,:,:);
q3 = qTraj(4,:,:);

dq0 = dq(1,:,:);
dq1 = dq(2,:,:);
dq2 = dq(3,:,:);
dq3 = dq(4,:,:);

wx = 2*(-q1.*dq0 + q0.*dq1 - q3.*dq2 + q2.*dq3);
wy = 2*(-q2.*dq0 + q3.*dq1 + q0.*dq2 - q1.*dq3);
wz = 2*(-q3.*dq0 - q2.*dq1 + q1.*dq2 + q0.*dq3);

w = [wx; wy; wz];

end

function wp = q2wp(qTraj, dt)

dq = derivative(qTraj, dt);

q0 = qTraj(1,:,:);
q1 = qTraj(2,:,:);
q2 = qTraj(3,:,:);
q3 = qTraj(4,:,:);

dq0 = dq(1,:,:);
dq1 = dq(2,:,:);
dq2 = dq(3,:,:);
dq3 = dq(4,:,:);

wxp = 2*(-q1.*dq0 + q0.*dq1 + q3.*dq2 - q2.*dq3);
wyp = 2*(-q2.*dq0 - q3.*dq1 + q0.*dq2 + q1.*dq3);
wzp = 2*(-q3.*dq0 + q2.*dq1 - q1.*dq2 + q0.*dq3);

wp = [wxp; wyp; wzp];

end

function msd = msd_fft(x)

if iscolumn(x)
  x = x.';
end

nComps = size(x, 1);
N = length(x);

D = zeros(nComps, N+1);
D(:,2:end) = x.^2;


% D = D.sum(axis=1)
% D = np.append(D,0)
S2 = runprivate('autocorrfft', x, 2, 0, 0, 0);

Q = 2*sum(D, 2);
S1 = zeros(nComps, N);

for m = 1:N
    Q = Q - D(:, m) - D(:, end-m);
    S1(:, m) = Q/((N+1)-m);
end

msd = S1 - 2*S2;

end