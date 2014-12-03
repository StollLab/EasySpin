function Gamma = rdogamma(basis,R,nSpin)


Lmax = max(basis(:,1));
L   = basis(:,1);
M   = basis(:,2);
K   = basis(:,3);
jK  = basis(:,4);
index = basis(:,5);
nSpace = sum((2*(0:Lmax)+1).^2);
nBasis = nSpace*nSpin;
nNonZero = length(basis);


%Gamma = sparse(zeros(nbasis,nbasis));


Rx = R(1,1);
Ry = R(2,2);
Rz = R(3,3);
Rd = 0.25*(Rx-Ry);
Rperp = 0.5*(Rx+Ry);

i = 1;
for iBasis = 1:nNonZero
  L1_  =  L(iBasis);
  M1_  =  M(iBasis);
  jK1_ = jK(iBasis);
  K1_  =  K(iBasis);
  if K1_ == 0
    deltaK1 = 1;
  else
    deltaK1 = 0;
  end
  idx1 = index(iBasis) - 1;
  
  for jBasis = iBasis:nNonZero
    L2_  =  L(jBasis);
    if (L1_ ~= L2_), break; end;
    L_ = L1_;
    M2_  =  M(jBasis);
    if (M1_ ~= M2_), continue; end;
    jK2_ = jK(jBasis);
    if (jK1_ ~= jK2_), continue; end;
    K2_  =  K(jBasis);
    if (abs(K1_-K2_) > 2), continue; end;
    if K2_ == 0
      deltaK2 = 1;
    else
      deltaK2 = 0;
    end
    idx2 = index(jBasis) - 1;
    
    if K1_ == K2_
      K_ = K1_;
      val_ = Rperp*(L1_*(L1_+1)-K_^2) + Rz*K_^2;
    else     
      if K2_ == K1_+2
        deltaUp = 1;
      else
        deltaUp = 0;
      end
      if K2_ == K1_-2
        deltaDown = 1;
      else
        deltaDown = 0;
      end
      
      Nup   = sqrt((L_-K2_-1)*(L_-K2_)*(L_+K2_+1)*(L_+K2_+2));
      Ndown = sqrt((L_+K2_-1)*(L_+K2_)*(L_-K2_+1)*(L_-K2_+2));
      NK = sqrt((1+deltaK1)*(1+deltaK2));
      val_ = Rd*NK*(Nup*deltaDown + Ndown*deltaUp); 
    end
      
      spinblock = eye(nSpin)*val_;
      [row,col,val] = find(spinblock);
      row = row + idx1;
      col = col + idx2;
      indices = i:i+numel(row)-1;
      
      bra(indices) = row;
      ket(indices) = col;
      el(indices)  = val;
      
      i = i + numel(row);
      
  end
end
    
Gamma = sparse(bra,ket,el,nBasis,nBasis);
Gamma = Gamma - (Gamma - diag(diag(Gamma))).';


return