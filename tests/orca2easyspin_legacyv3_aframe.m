function ok = test()

% Legacy ORCA 3.0.3 files: hyperfine frames of a larger molecule (Pr_EPR)
% from the binary property file, compared with reference values, and
% agreement of hyperfine tensors between property file and main output file

folder = 'orca/v3.0.3/';

T = @(v,ang) erot(ang).'*diag(v)*erot(ang);
maxdiff = @(X,Y) max(abs(X(:)-Y(:)));

SysP = orca2easyspin([folder 'Pr_EPR.prop']);
AFrame_ref = [
   2.704188721737625   0.342303276412723   3.804135252240579
   5.703961547204178   1.758814960534842   1.285881496415190
   4.367658953893748   1.016119233116642   0.924490911561192
   5.224528319045627   1.119158243324929   1.034825024969174
   4.946659501467618   0.281003774604616   5.501265355009233
   4.250629412404098   0.246486610983774   0.990357494553786
   1.949808589552022   1.923345926558176   3.007831332340492
   ];
ok(1) = areequal(SysP.AFrame,AFrame_ref,1e-10,'abs');

SysM = orca2easyspin([folder 'Pr_EPR.out']);
okk = strcmp(SysM.Nucs,SysP.Nucs) && isequal(SysM.NucsIdx,SysP.NucsIdx);
for n = 1:numel(SysP.NucsIdx)
  okk = okk && maxdiff(T(SysM.A(n,:),SysM.AFrame(n,:)),T(SysP.A(n,:),SysP.AFrame(n,:)))<0.01;  % MHz
end
ok(2) = okk;
