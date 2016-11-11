function [B] = PoissonRetinoMap(Response,ScreenDim)
%PoissonRetinoMap.m
%   Take data from an LFP retinotopic mapping experiment and find the
%    retinotopic center of mass. 
%INPUT: Response - matrix with stimulus locations and an indicator function
%        for hit (1) or miss (0) in response to each stimulus
%       ScreenDim - dimensions of the screen being used, in pixels,
%        as [horizontal pixels,vertical pixels]
%OUTPUT: Beta - a number of parameters-by-1 vector that
%         contains the fit parameter values for the model
%
%Created: 2016/11/10, 24 Cummington Mall
% Byron Price
%Updated: 2016/11/10
%  By: Byron Price

Y = squeeze(Response(:,3));

DesignMat = [Response(:,1) Response(:,2) Response(:,1).*Response(:,1) ...
    Response(:,2).*Response(:,2) Response(:,1).*Response(:,2)];

[B,dev,stats] = glmfit(DesignMat,Y,'poisson');

Standard_Error = stats.se;
T = table(B,Standard_Error,'RowNames',{'Intercept','x','y','x^2','y^2','xy'});
display(T);
display(sprintf('Model Deviance: %3.2f',dev));

xPos = 1:ScreenDim(1);yPos = 1:ScreenDim(2);
[xPos,yPos] = meshgrid(xPos,yPos);
modelFit = exp(B(1)+B(2).*xPos+B(3).*yPos+B(4).*xPos.*xPos+B(5).*yPos.*yPos+B(6).*xPos.*yPos);

figure();imagesc(modelFit);colormap('jet');colorbar;
title('Retinotopic Map Fit for Poisson GLM');
end

