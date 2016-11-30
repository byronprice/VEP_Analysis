function [] = LikelyRetinoMap(X,VEPs,VEP_Prototype,ScreenDim)
%LikelyRetinoMap.m
%   Take data from an LFP retinotopic mapping experiment and find the
%    retinotopic center of mass. 
%INPUT: X - a number of stimuli-by-2 matrix of the (x,y) positions of the
%        stimuli presented to the animal
%       VEPs - a number of stimuli-by-200 matrix of VEPs recorded in
%        response to the presented stimuli, 200 data points sampled at 
%        1000 Hz
%       ScreenDim - dimensions of the screen being used, in pixels,
%        as [horizontal pixels,vertical pixels]
%OUTPUT: B - a number of parameters-by-1 vector that
%         contains the fit parameter values for the model
%
%Created: 2016/11/11, 24 Cummington Mall
% Byron Price
%Updated: 2016/11/11 
%  By: Byron Price

numVEPs = length(X);
vepLen = length(VEP_Prototype);
DistFun = @(stimCenter,centerVals) (ceil(sqrt((stimCenter(1)-centerVals(:,1)).^2+(stimCenter(2)-centerVals(:,2)).^2))+1);

xPos = 1:ScreenDim(1);yPos = 1:ScreenDim(2);
lenx = length(xPos);leny = length(yPos);
centerVals = zeros(lenx*leny,2);

count = 1;
for ii=1:lenx
    for jj=1:leny
        centerVals(count,1) = xPos(ii);
        centerVals(count,2) = yPos(jj);
        count = count+1;
    end
end

count = 1;
Prototype = repmat(VEP_Prototype,[numVEPs,1]);
Y = zeros(numVEPs*vepLen,1);
bigX = zeros(numVEPs*vepLen,2);
for ii=1:numVEPs   
    Y(count:count+vepLen-1) = VEPs(ii,:)';
    bigX(count:count+vepLen-1,:) = repmat(X(ii,1:2),[vepLen,1]);
    count = count+vepLen;
end

% calculate likelihood (p(VEP | distance from center of mass) = normal)
b = [1,250,0];
maxDist = ceil(sqrt((xPos(end)-xPos(1)).^2+(yPos(end)-yPos(1)).^2));
likelyFun = @(b,distToCenterMass) (b(1)*exp(-(distToCenterMass.^2)./(2*b(2)*b(2)))+b(3));
GaussianModel = likelyFun(b,(0:maxDist)');

sigmasquare = 500;
loglikely = zeros(lenx*leny,length(sigmasquare));
for ii=1:numVEPs*vepLen
    stimCenter = [bigX(ii,1),bigX(ii,2)];
    allPossDists = DistFun(stimCenter,centerVals);
    RHS = GaussianModel(allPossDists);
    squaredDiff = (Y(ii)-RHS.*Prototype(ii)).^2;
    for jj=1:length(sigmasquare)
        loglikely(:,jj) = loglikely(:,jj)+log((1./sqrt(2*pi*sigmasquare(jj))).*exp(-squaredDiff./(2*sigmasquare(jj))));
    end
end

[maxVals,rowInds] = max(loglikely);
[~,columnInd] = max(maxVals);
rowInd = rowInds(columnInd);

centerX_hat = centerVals(rowInd,1);
centerY_hat = centerVals(rowInd,2);
sigmasquare_hat = sigmasquare(columnInd);

loglikely = loglikely(:,columnInd);
loglikely = reshape(loglikely,[leny,lenx]);
figure();subplot(3,1,1);imagesc(xaxis,yaxis,loglikely);set(gca,'YDir','normal');
colorbar;
end

%  gradient ascent for max likelihood
% x=500;y=500;sigmasquare = 100;
% eps = [-0.5,0,0.5];
% theta = [x,y,sigmasquare];
% for ii=1:numIterations
%    loglikely = zeros(length(theta),3);
%    for jj=1:numVEPs*vepLen
%        for kk=1:length(theta)
%            tempTheta = theta;
%            for ll=1:3
%                tempTheta(kk) = tempTheta(kk)+eps(ll);
%                stimCenter = [bigX(ii,1),bigX(ii,2)];
%                dist = DistFun(stimCenter,tempTheta(1:2));
%                RHS = GaussianModel(dist);
%                squaredDiff = (Y(ii)-RHS.*Prototype(ii)).^2;
%                loglikely(kk,ll) = loglikely(kk,ll)+log((1./sqrt(2*pi*tempTheta(3))).*exp(-squaredDiff./(2*tempTheta(3))));
%            end
%        end
%    end
%    
% end