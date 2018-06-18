% cardamom_rejectionsample3d  Draw random samples from a 3-variate
%                             distribution using the rejection method.
%
%   [xSamples, ySamples, zSamples] = ...
%            cardamom_rejectionsample3d(pdf, xGrid, yGrid, zGrid, nSamples)
%
%   Input:
%     pdf            numeric array, size = (N,M,P)
%                    probability distribution function
%
%     xGrid          numeric array, size = (1,N)
%                    grid of x-values in the domain of the pdf
%
%     yGrid          numeric array, size = (1,N)
%                    grid of y-values in the domain of the pdf
%
%     zGrid          numeric array, size = (1,N)
%                    grid of z-values in the domain of the pdf
%
%     nSamples       integer
%                    number of samples to draw from the pdf
%
%
%   Output:
%     xSamples       numeric array, size = (1,nSamples)
%                    samples drawn from the x-domain of the pdf
%
%     ySamples       numeric array, size = (1,nSamples)
%                    samples drawn from the y-domain of the pdf
%
%     zSamples       numeric array, size = (1,nSamples)
%                    samples drawn from the z-domain of the pdf

function [xSamples, ySamples, zSamples] = cardamom_rejectionsample3d(pdf, xGrid, yGrid, zGrid, nSamples)

NxBins = size(xGrid);
NyBins = size(yGrid);
NzBins = size(zGrid);

c = max(pdf(:));

iSample = 1;
xSamples = zeros(1,nSamples);
ySamples = zeros(1,nSamples);
zSamples = zeros(1,nSamples);
while iSample < nSamples + 1
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
end

end