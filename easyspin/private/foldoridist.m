% Fold orientational distribution function into unorifunique region of grid
% symmetry.
%
% Input:
%   f        ... orientational distribution function, f(beta,gamma)
%                where (alpha,beta,gamma) are the Euler angles for the
%                transformation from sample frame to molecular frame
%   GridSymm ... grid symmetry group ('D2h', 'C2h', 'Ci', 'C1')
%
% Output:
%   orifun   ... folded orientational distribution function, orifun(beta,gamma)

function orifun = foldoridist(f,GridSymmetry)
switch GridSymmetry
  case 'O3'
    ff = @(beta,gamma)f(beta,gamma).*sin(beta);
    val = integral2(ff,0,pi,0,2*pi); % integral over sphere
    orifun = @(beta,gamma) val;
  case 'Dinfh'
    orifun = @(beta,gamma) gammaint(f,beta); % integral over all gamma
  case 'D2h'
    orifun = @(beta,gamma)...
      f(beta,gamma)+...       % E ( = identity)
      f(beta,pi-gamma)+...    % sigma(yz)
      f(beta,pi+gamma)+...    % C2(z)
      f(beta,2*pi-gamma)+...  % sigma(xz)
      f(pi-beta,gamma)+...    % sigma(xy)
      f(pi-beta,pi-gamma)+... % C2(y)
      f(pi-beta,pi+gamma)+... % i
      f(pi-beta,2*pi-gamma);  % C2(x)
  case 'C2h'
    orifun = @(beta,gamma)...
      f(beta,gamma)+...       % E
      f(beta,pi+gamma)+...    % C2(z)
      f(pi-beta,gamma)+...    % sigma(xy)
      f(pi-beta,pi+gamma);    % i
  case 'Ci'
    orifun = @(beta,gamma)...
      f(beta,gamma)+...       % E
      f(pi-beta,gamma+pi);    % i
  case 'C1'
    orifun = @(beta,gamma)...
      f(beta,gamma);          % E
  otherwise
    error('Orientational distribution folding are not supported for grid symmetry ''%s''.',GridSymmetry);
end
end

function v = gammaint(f,beta)
for k = numel(beta):-1:1
  v(k) = integral(@(gamma)f(beta(k),gamma)+f(pi-beta(k),gamma),0,2*pi); 
end
end
