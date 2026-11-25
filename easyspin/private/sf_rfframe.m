% sf_rfframe  Nuclear RF frame and rotating-wave error for RF pulses
%
%   [rfFrame,rwaError,Xp] = sf_rfframe(H,V,X0,Spins,nElectrons,drive)
%
%   RF pulses are propagated in a frame that rotates about the generator
%   G = V*diag(f)*V', which is diagonal in the eigenbasis V of H, so that H
%   commutes with G exactly and is not truncated (see rfFramePropagators in
%   sf_propagate). Each eigenstate k (column k of V, assigned to product
%   state k by sf_eigenbasis) gets the label f(k) = sum over nuclei i of
%   s_i*m_i(k), where m_i(k) is the m of nucleus i in that product state,
%   and s_i = +1 or -1 is the sense of precession of nucleus i in the
%   electron mS manifold of k: the sign of tr(P_mS*Iz_i*H), the first-order
%   z frequency of nucleus i in the manifold (+1 if zero). The sense can
%   differ between the manifolds, e.g. if the hyperfine coupling is larger
%   than twice the nuclear Zeeman frequency. With these labels, the energy
%   of the states rises with f along the NMR transitions, so all of them
%   rotate with the frame.
%
%   Rotating-wave approximation: the RF field 2*Re(c*exp(2i*pi*F*t))*X0 is
%   kept only where it rotates with the frame, i.e. the term c*X(k,l) where
%   f(k)-f(l) = -1, and its Hermitian conjugate (X = V'*X0*V). Relative to
%   H, a dropped term c*X(k,l) oscillates at E(k)-E(l)+F. Its leading
%   effect is a Bloch-Siegert shift, which, relative to the nutation
%   frequency, is the ratio of its amplitude |c*X(k,l)| to the frequency
%   distance |E(k)-E(l)+F|. The error of a pulse is the largest ratio over
%   all dropped terms, with F anywhere in the frequency band of the pulse.
%
%   RF field directions: X0 holds the RF operators of one or two field
%   directions. With two, n1 and n2, the RF field can point along any
%   direction u1*n1 + u2*n2 (u real, u1^2+u2^2 = 1), with the operator
%   u1*X0{1} + u2*X0{2}. The error is the largest over all these
%   directions, so that the choice of stepped pulses is the same for all of
%   them (a smooth signal for the average over the RF field direction).
%   For the elements a and b of a dropped term in X{1} and X{2}, the
%   largest amplitude |u1*a + u2*b| is the square root of the largest
%   eigenvalue of real([a b]'*[a b]). Xp is linear in X0, so Xp of a
%   direction is u1*Xp{1} + u2*Xp{2}.
%
%   Input:
%     H           Hamiltonian (MHz) in the simulation frame, commuting with
%                 the total electron Sz
%     V           eigenvectors of H, from sf_eigenbasis
%     X0          RF excitation operators (MHz/mT), product basis, cell
%                 array with one or two operators: the field B2 (mT) along
%                 the direction of X0{k} gives B2*X0{k}
%     Spins       spin quantum numbers (Sys.Spins), electrons first
%     nElectrons  number of electron spins
%     drive       RF pulses, one row per pulse and point: largest field
%                 amplitude (mT), lowest and highest frequency of the band
%                 (MHz); zeros(0,3) if there are none
%
%   Output:
%     rfFrame     structure with the fields
%       V           eigenvectors of H (columns), as input
%       E           eigenvalues (MHz), column
%       f           labels, column
%                 The caller adds Xp for the RF field direction (see
%                 rfFramePropagators in sf_propagate).
%     rwaError    relative error of the rotating-wave approximation for each
%                 row of drive, the largest over the RF field directions,
%                 column (Inf if a dropped term is resonant)
%     Xp          for each operator in X0, the elements of X with
%                 f(k)-f(l) = +1, zero elsewhere, cell array

function [rfFrame,rwaError,Xp] = sf_rfframe(H,V,X0,Spins,nElectrons,drive)

E = real(diag(V'*H*V));

% Electron manifold of each product state
mS = 0;
for iSpin = 1:nElectrons
  mS = mS + full(diag(sop(Spins,[iSpin,3])));
end
[~,~,manifold] = unique(mS);

% Labels: m of each nucleus, signed by its sense of precession in the
% manifold, tr(P_mS*Iz*H) = sum of m.*diag(H) over the manifold
Hdiag = real(diag(H));
f = zeros(size(H,1),1);
for iSpin = nElectrons+1:numel(Spins)
  m = full(diag(sop(Spins,[iSpin,3])));
  for iManifold = 1:max(manifold)
    inManifold = manifold==iManifold;
    s = sign(sum(m(inManifold).*Hdiag(inManifold)));
    if s==0, s = 1; end
    f(inManifold) = f(inManifold) + s*m(inManifold);
  end
end

rfFrame.V = V;
rfFrame.E = E;
rfFrame.f = f;
kept = abs(f-f.'-1)<1e-9;
X = cell(size(X0));
Xp = cell(size(X0));
for k = 1:numel(X0)
  X{k} = V'*X0{k}*V;
  Xp{k} = X{k}.*kept;
end

% Largest amplitude of each element over the RF field directions
if numel(X)==1
  Xabs = abs(X{1});
else
  p = abs(X{1}).^2;
  q = abs(X{2}).^2;
  r = real(conj(X{1}).*X{2});
  Xabs = sqrt((p+q)/2 + sqrt(((p-q)/2).^2 + r.^2));
end

% Rotating-wave error of each pulse: ratio of the amplitude of each dropped
% term to its frequency distance from the band of the pulse
dropped = abs(f-f.'+1)>1e-9;
maxAmp = max(Xabs(:));
Xabs = Xabs.*dropped;
Xabs(Xabs<=1e-12*maxAmp) = 0;
resonance = E.' - E;  % E(l)-E(k): frequency at which the term c*X(k,l) is static
rwaError = zeros(size(drive,1),1);
for r = 1:size(drive,1)
  distance = max(max(drive(r,2)-resonance,resonance-drive(r,3)),0);
  ratio = drive(r,1)*Xabs./distance;
  ratio(Xabs==0) = 0;
  rwaError(r) = max([0; ratio(:)]);
end

end
