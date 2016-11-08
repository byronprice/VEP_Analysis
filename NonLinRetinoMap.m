 function [MapFun,Beta] = NonLinRetinoMap(X,VEPs,ScreenDim)
%NonLinRetinoMap.m
%   Take data from an LFP retinotopic mapping experiment and find the
%    retinotopic center of mass. 
%INPUT: X - a number of stimuli-by-2 matrix of the (x,y) positions of the
%        stimuli presented to the animal
%       VEPs - a number of stimuli-by-200 matrix of VEPs recorded in
%        response to the presented stimuli, 200 data points sampled at 
%        1000 Hz
%       ScreenDim - dimensions of the screen being used, in pixels,
%        as [horizontal pixels,vertical pixels]
%OUTPUT: MapFun - function that represents the model of the retinotopic
%         map, i.e. a 2D Gaussian
%        Beta - a number of parameters-by-1 vector that
%         contains the fit parameter values for the model
%
%Created: 2016/10/18, 725 Commonwealth Avenue
% Byron Price
%Updated: 2016/10/18
%  By: Byron Price

load('VEP_Prototype.mat');

numTests = 50;
b0 = [random('Normal',ScreenDim(1)/2,ScreenDim(1)/5,[numTests,1]),...
    random('Normal',ScreenDim(2)/2,ScreenDim(2)/5,[numTests,1]),random('Normal',50000,5000,[numTests,1]),...
    random('Normal',20,4,[numTests,1]),random('Normal',50000,5000,[numTests,1]),...
    random('Normal',1,0.2,[numTests,1])];

Beta = zeros(numTests,6);
for ii=1:numTests
    tempX = X;
    tempY = VEPs;
    for jj=1:10
        index = random('Discrete Uniform',length(tempX),1);
        tempX(index,:) = [];
        tempY(index,:) = [];
    end

    lb = [0,0,0,0,0,0];
    ub = [ScreenDim(1),ScreenDim(2),1e7,1e5,1e7,1e2];
    vepLen = length(VEP_Prototype);
    bigVEPs = repmat(VEP_Prototype',[length(tempX),1]);
    myFun = @(b) bigVEPs.*repmat(b(6).*exp(-0.5.*((b(5)/(b(3)*b(5)-b(4)*b(4))).*(tempX(:,1)-b(1)).^2-...
        (tempX(:,1)-b(1)).*(tempX(:,2)-b(2)).*(2*b(4)/(b(3)*b(5)-b(4)*b(4)))+...
        (b(3)/(b(3)*b(5)-b(4)*b(4))).*(tempX(:,2)-b(2)).^2)),[1,vepLen])-tempY;
    
    [b,~] = lsqnonlin(myFun,b0(ii,:),lb,ub);
    Beta(ii,:) = b;
end
MapFun = @(b,x) b(6).*exp(-0.5.*((b(5)/(b(3)*b(5)-b(4)*b(4))).*(x(:,1)-b(1)).^2-...
             (x(:,1)-b(1)).*(x(:,2)-b(2)).*(2*b(4)/(b(3)*b(5)-b(4)*b(4)))+...
                   (b(3)/(b(3)*b(5)-b(4)*b(4))).*(x(:,2)-b(2)).^2));

trueBeta = median(Beta,1);
display(sprintf('Center at [X,Y]: %3.3f %3.3f',trueBeta(1),trueBeta(2)));
x = 0:10:2500;
y = 0:10:1500;
Z = zeros(length(x),length(y));
for ww=1:length(x)
   for xx=1:length(y)
       Z(ww,xx) = MapFun(trueBeta,[x(ww),y(xx)]);
   end
end
figure();imagesc(x,y,Z');set(gca,'YDir','normal');

Beta

% DataPoints = length(VEPs);
% vepLen = length(VEP_Prototype);
% Response = zeros(DataPoints*vepLen,1);
% Predictors = zeros(DataPoints*vepLen,2);
% 
% index = 1;
% for ii=1:DataPoints
%    Response(index:index+vepLen-1,1) = VEPs(ii,:)'./VEP_Prototype;
%    Predictors(index:index+vepLen-1,:) = repmat(X(ii,:),[vepLen,1]);
%    index = index+vepLen;
% end
% 
% model = fitnlm(Predictors,Response,MapFun,b0(end,:),'CoefficientNames',{'CenterX','CenterY','SigmaX','SigmaXY','SigmaY','ScaleFactor'})
end

