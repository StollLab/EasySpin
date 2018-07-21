function [Gamma] = s_relaxationsuperoperator(System)
%  [gamma] = changeup(system)
%
%  The function changeup calculates the relaxation superoperators for a 
%  given spinsystem system.sqn, which is a vector of spin quantum numbers.
%  Relaxation times can either be given as a single or in matrix form. 
%  In the former case, the same value is applied all T1 or T2 relaxation
%  pathways:
%  
%  System.T1 = 1000;
%  System.T2 = 500;
%
%  The latter allows to define T1 and T2 values to all transitions 
%  seperately. The used matrix has the same structure as a denisty matrix
%  of the system, with elements corresponding to the same states as in
%  sigma. The off-diagonal elements correlate the pathways between the
%  states in the first dimension and the second dimension. Therefore it is
%  necessary to define the upper triangle only.
% 
%  Assign different relaxation times to all transitions in above example:
% 
%                aa  ab     ba      bb
%  System.T1 =  [0   10^7   10^3    10^4;  aa
%                0   0      10^4    10^3;  ab
%                0   0      0       10^7;  ba
%                0   0      0       0   ;] bb
% 
%                aa  ab     ba      bb
%  System.T2 =  [0   10^3   10^2    10^3;  aa
%                0   0      10^3    10^2;  ab
%                0   0      0       10^5;  ba
%                0   0      0       0   ;] bb
% 
%  From this, the program builds the relaxation superator gamma in 
%  Liouville space.
%
%  St. Pribitzer, 2017

[t11, t12] = size(System.T1);
[t21, t22] = size(System.T2);

HilbertDimension = prod(System.Spins*2+1);

T1ud = 0;

% checks the format of the input and expands it to matrices if necessary
% for T1 times, which allows for easier processing
if t11 == 1 && t12 == 1
    T1 = System.T1*ones(HilbertDimension); 
else
    T1 = System.T1;
    if any(any(tril(T1)))
       T1ud = 1; 
    end
end

% checks the format of the input and expands it to matrices if necessary
% for T2 times, which allows for easier processing
if t21 ~= 1 || t22 ~= 1
    T2 = System.T2;
%     if ~any(any(tril(T2)))
        % mirrors upper triangle of matrix to lower triangle
        T2 = triu(T2) + triu(T2,1)';
%     end
elseif t21 == 1 && t22 == 1 
    T2 = System.T2*ones(HilbertDimension);   
end



% Assigns default values or inf to all T1 relaxation paths which are zero.
T1(T1==0) = 1e10;
T2(T2==0) = 1e10;

Gamma = zeros((HilbertDimension)^2);

kk = 1;
jj = 2;

% calculates the positions of longitudinal relaxation, populations are
% never on the diagonal elements of the relaxation superoperator.
for k1 = 1:HilbertDimension+1:HilbertDimension^2
    
    for k2 = k1+HilbertDimension+1 : HilbertDimension+1 : HilbertDimension^2
        
        Gamma(k1,k2) = -1/T1(kk,jj);
        
        if T1ud
            Gamma(k2,k1) = -1/T1(jj,kk);
        else
            Gamma(k2,k1) = -1/T1(kk,jj);
        end
        % since we are using only reduced density matrices (difference of
        % the equilibrium state to the current state), it is not necessary
        % to compensate the T1 decay of off-diagonal elements/states with
        % the diagonal elements:
%         gamma(k1,k1)=gamma(k1,k1)+1/T1(kk,jj);
%         gamma(k2,k2)=gamma(k2,k2)+1/T1(kk,jj);
        jj = jj+1;
        
    end
    
    kk = kk+1;
    jj = kk+1;
    
end

% populates positions of transverse relaxation elements
n = 1;
% rewrites T2 matrix into a vector, to allow effective processing
T2vec = reshape(T2,HilbertDimension*HilbertDimension,1);

for k = 2:HilbertDimension^2-1
    
    if k ~= n+1+n*HilbertDimension
        Gamma(k,k) = 1/T2vec(k);
    else
        n = n+1;
    end 
end

end

