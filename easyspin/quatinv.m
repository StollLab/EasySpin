%  quatinv Invert a quaternion.
%
%  qinv = quatvecmult(q);
%
%  Input:
%     q          4x... array, normalized quaternion
%
%  Output:
%     qinv       4x... array, normalized quaternion

function qinv = quatvecmult(q)
    
if size(q, 1) ~= 4
  error('Size of first dimension of the quaternion must equal 4.')
end

if any(1.0-sum(q.*q, 1) > 1e-5)
    error('Input quaternion is not normalized.')
end

inverter = [1; -1; -1; -1];

qinv = bsxfun(@times, inverter, q);

end