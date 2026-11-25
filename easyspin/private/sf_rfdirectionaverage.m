% sf_rfdirectionaverage  Average of a signal over the RF field direction
%
%   [avg,nDirections,converged] = sf_rfdirectionaverage(signal,chiRange)
%
%   The RF coil field is linearly polarized, along y(Lab) in crystals
%   (perpendicular to the microwave field along x(Lab)). Powder grids have
%   only two angles (phi,theta), with the third Euler angle chi = 0.
%   Rotating a molecule by chi about the static field is equivalent to
%   rotating the RF field direction by -chi with the molecule fixed: the
%   microwave pulses only change phase, which cancels against the phase of
%   the detection, since the secular Hamiltonian, the drift part of the
%   microwave pulses and the RF operator commute with the total electron
%   Sz. This holds in a shifted simulation frame. In the lab frame, it is
%   approximate, as is the chi = 0 grid itself. So in powders the signal is
%   averaged over the RF field direction n = cos(chi)*x(Lab) +
%   sin(chi)*y(Lab), with chi in [0,chiRange): pi with a random RF phase (n
%   and -n differ only by the RF phase pi), 2*pi with a phase-locked RF.
%
%   The integrand is smooth and periodic, so the trapezoid rule converges
%   exponentially. Starting with 4 directions, the number is doubled until
%   the averages over N directions and over its subset of N/2 differ by
%   less than 1e-5 of the largest signal, up to 256 directions. For weak
%   RF pulses, the signal is linear in the squared RF matrix elements,
%   which have harmonics up to 2*chi, so 4 directions suffice there.
%
%   The average is needed also without the hyperfine enhancement: the RF
%   matrix elements depend on the direction through the tilt of the
%   nuclear quantization axes from the static field, which is large near
%   the cancellation condition A = 2*nu_I (errors of tens of percent of
%   the ENDOR effect with the RF field along y(Lab) only).
%
%   Input:
%     signal       function handle, signal(chi) returns the signal (numeric
%                  array) for the RF field direction chi (radians from
%                  x(Lab)); called in the order of the directions, starting
%                  with chi = 0
%     chiRange     range of chi (radians): pi or 2*pi
%
%   Output:
%     avg          signal averaged over the directions
%     nDirections  number of directions used
%     converged    false if the average did not converge

function [avg,nDirections,converged] = sf_rfdirectionaverage(signal,chiRange)

tolerance = 1e-5;
maxDirections = 256;

% sums over all directions (sumAll) and over every second one (sumHalf)
nDirections = 4;
chi = (0:nDirections-1)*chiRange/nDirections;
sumAll = 0;
sumHalf = 0;
for j = 1:nDirections
  s = signal(chi(j));
  sumAll = sumAll + s;
  if mod(j,2)==1, sumHalf = sumHalf + s; end
end

converged = true;
while true
  avg = sumAll/nDirections;
  change = max(abs(avg(:) - sumHalf(:)/(nDirections/2)));
  if change<=tolerance*max(abs(avg(:))), break; end
  if nDirections>=maxDirections
    converged = false;
    break
  end
  % add the directions halfway between the earlier ones
  sumHalf = sumAll;
  for chiNew = chi + chiRange/(2*nDirections)
    sumAll = sumAll + signal(chiNew);
  end
  nDirections = 2*nDirections;
  chi = (0:nDirections-1)*chiRange/nDirections;
end

end
