function [xSamples, ySamples, zSamples] = rejectionsample3d(pdf, xGrid, yGrid, zGrid, nSamples)

NxBins = size(xGrid);
NyBins = size(yGrid);
NzBins = size(zGrid);

c = max(max(max(pdf)));

iSample = 1;
count = 0;
xSamples = zeros(1,nSamples);
ySamples = zeros(1,nSamples);
zSamples = zeros(1,nSamples);
while iSample < nSamples
  xProposal = randi(NxBins);
  yProposal = randi(NyBins);
  zProposal = randi(NzBins);
  q = c;
  if rand() < pdf(xProposal,yProposal,zProposal)/q
    xSamples(iSample) = xGrid(xProposal);
    ySamples(iSample) = yGrid(yProposal);
    zSamples(iSample) = zGrid(zProposal);
    iSample = iSample + 1;
  end
  count = count + 1;
end

end