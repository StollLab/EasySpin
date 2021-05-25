% Fold orientational distribution function into unique region of grid
% symmetry.
%
% Input:
%   f        ... orientational distribution function, f(phi,theta)
%   GridSymm ... grid symmetry group ('D2h', 'C2h', 'Ci', 'C1')
%
% Output:
%   orifun   ... folded orientational distribution function, orifun(phi,theta)

function orifun = foldoridist(f,GridSymmetry)
switch GridSymmetry
  case 'O3'
    ff = @(phi,theta)f(phi,theta).*sin(theta);
    val = integral2(ff,0,2*pi,0,pi); % integral over sphere
    orifun = @(phi,theta) val;
  case 'Dinfh'
    orifun = @(phi,theta) phiint(f,theta); % integral over all phi
  case 'D2h'
    orifun = @(phi,theta)...
      f(phi,theta)+...       % E ( = identity)
      f(pi-phi,theta)+...    % sigma(yz)
      f(pi+phi,theta)+...    % C2(z)
      f(2*pi-phi,theta)+...  % sigma(xz)
      f(phi,pi-theta)+...    % sigma(xy)
      f(pi-phi,pi-theta)+... % C2(y)
      f(pi+phi,pi-theta)+... % i
      f(2*pi-phi,pi-theta);  % C2(x)
  case 'C2h'
    orifun = @(phi,theta)...
      f(phi,theta)+...       % E
      f(pi+phi,theta)+...    % C2(z)
      f(phi,pi-theta)+...    % sigma(xy)
      f(pi+phi,pi-theta);    % i
  case 'Ci'
    orifun = @(phi,theta)...
      f(phi,theta)+...       % E
      f(phi+pi,pi-theta);    % i
  case 'C1'
    orifun = ExpOrdering;    % E
  otherwise
    error('Orientational distribution folding are not supported for grid symmetry ''%s''.',GridSymmetry);
end
end

function v = phiint(f,theta)
for k = numel(theta):-1:1
  v(k) = integral(@(phi)f(phi,theta(k))+f(phi,pi-theta(k)),0,2*pi); 
end
end
