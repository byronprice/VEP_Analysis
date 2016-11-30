function [B] = GLMRetinoMap(X,VEPs,VEP_Prototype,ScreenDim)
%GLMRetinoMap.m
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
%Created: 2016/11/10, 24 Cummington Mall
% Byron Price
%Updated: 2016/11/10
%  By: Byron Price

numVEPs = length(X);
vepLen = length(VEP_Prototype);

numCoeffs = 5;
DesignMat = zeros(numVEPs*vepLen,numCoeffs);
Y = zeros(numVEPs*vepLen,1);

% create design matrix
count = 1;
for ii=1:numVEPs
    temp = repmat(X(ii,:),[vepLen,1]);
    for jj=1:numCoeffs
        DesignMat(count:count+vepLen-1,jj) = temp(:,jj).*VEP_Prototype;
    end
    Y(count:count+vepLen-1) = VEPs(ii,:)';
    count = count+vepLen;
end

[B,dev,stats] = glmfit(DesignMat,Y,'normal');

Standard_Error = stats.se;
P_Values = stats.p;
T = table(B,Standard_Error,P_Values,'RowNames',{'Intercept','x','y','xx','yy','xy'});
display(T);
display(sprintf('Model Deviance: %3.2f',dev));

vepProtoMin = abs(min(VEP_Prototype));
xPos = 1:ScreenDim(1);yPos = 1:ScreenDim(2);
[xPos,yPos] = meshgrid(xPos,yPos);
modelFit = B(1)+vepProtoMin.*(B(2).*xPos+B(3).*yPos+B(4).*xPos.*xPos+B(5).*yPos.*yPos+B(6).*xPos.*yPos);
figure();imagesc(modelFit);h = colorbar;ylabel(h,'VEP Negativity');
title('Retinotopic Map Fit for Normal GLM');
end

