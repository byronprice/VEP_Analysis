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

Y = Response(:,3);

% xcen = mean(Response(:,1));
% ycen = mean(Response(:,1));
% Response(:,1) = Response(:,1)-xcen;
% Response(:,2) = Response(:,2)-ycen;

Response(:,1) = Response(:,1)./max(Response(:,1));
Response(:,2) = Response(:,2)./max(Response(:,2));
DesignMat = [Response(:,1) Response(:,2) Response(:,1).*Response(:,1) ...
    Response(:,2).*Response(:,2) Response(:,1).*Response(:,2)];

%N = ones(length(DesignMat),1);
%[B,dev,stats] = glmfit(DesignMat,[Y,N],'binomial');
[B,dev,stats] = glmfit(DesignMat,Y,'poisson');

Standard_Error = stats.se;
P_Values = stats.p;
T = table(B,Standard_Error,P_Values,'RowNames',{'Intercept','x','y','xx','yy','xy'});
display(T);
display(sprintf('Model Deviance: %3.2f',dev));


xaxis = 0:0.01:1;yaxis = 0:0.01:1;
totalPos = length(xaxis)*length(yaxis);
X = zeros(totalPos,5);

count = 1;
for ii=1:length(xaxis)
    for jj=1:length(yaxis)
        X(count,:) = [xaxis(ii)' yaxis(jj)' xaxis(ii)'.*xaxis(ii)' yaxis(jj)'.*...
            yaxis(jj)' xaxis(ii)'.*yaxis(jj)'];
        count = count+1;
    end
end

modelFit = glmval(B,X,'log');

%modelFit = exp(B(1)+B(2).*xPos+B(3).*yPos+B(4).*xPos.*xPos+B(5).*yPos.*yPos+B(6).*xPos.*yPos)./(1+...
  %  exp(B(1)+B(2).*xPos+B(3).*yPos+B(4).*xPos.*xPos+B(5).*yPos.*yPos+B(6).*xPos.*yPos));

modelFit = reshape(modelFit,[length(yaxis),length(xaxis)]);
figure();imagesc(xaxis,yaxis,modelFit);colorbar;set(gca,'YDir','normal');
title('Retinotopic Map Fit for Poisson GLM');
end

