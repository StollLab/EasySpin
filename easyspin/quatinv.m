%  quatinv Invert a unit quaternion.
%
%  qinv = quatvecmult(q);
%
%  Input:
%     q              numeric, size = (4,...)
%                    normalized quaternion
%
%  Output:
%     qinv           numeric, size = (4,...)
%                    inverse of normalized quaternion

function qinv = quatinv(q)
    
if size(q, 1) ~= 4 || ~isnumeric(q)
  error('Input must be numeric with size (4,...).')
end

inverter = [1; -1; -1; -1];

qinv = bsxfun(@times, inverter, q);

end