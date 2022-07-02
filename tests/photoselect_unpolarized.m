function ok = test()

% Test the photoselection weight calculated for an unpolarized beam is consistent
% with an explicit average of photoselection weights over E-field directions of
% polarized beams.

tdm = [32 134]*pi/180;
ori = [77 11 222]*pi/180;
k = 'y';

% Reference value calculated by photoselect()
w_ref = photoselect(tdm,ori,k,NaN);

% Approximate explicit integral over polarization angle alpha
nAlpha = 1000;
alpha = linspace(0,pi,nAlpha);
w_test = zeros(1,nAlpha);
for i = 1:nAlpha
  w_test(i) = photoselect(tdm,ori,k,alpha(i));
end
w_test = trapz(w_test/nAlpha);

ok = areequal(w_test,w_ref,1e-3,'rel');

end
