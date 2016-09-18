% VEP_Likelihood.m
%  Get the likelihood function L(THETA | data), where theta is the distance
%   to the retinotopic center of mass and the data is the probability of a 
%   VEP-like event at that distance ... so, at the center of mass, maybe
%   the probability is 0.7, but much further away, it's 0.2 . The VEP-like
%   event is defined as a peak negativity < -150 with latency to peak
%   negativity between 60 and 120 msec on an individual trial. The idea is
%   that, under the hypothesis that the current position is the center of
%   mass, I would expect a probability of 0.7, so if I find such a
%   probability, I'm likely to be at the center of mass. If I'm far away
%   from the center of mass, then I expect a probability of 0.2 .

fileList = dir('RetinoMap*_*.mat');
numFiles = length(fileList);

fullWin = 1:200;
minWin = [60,120]; % milliseconds
maxWin = 1:250;
Pr = zeros(3000,2);
totResponse = cell(numFiles,1);
totnumChans = zeros(numFiles,1);
totnumStimuli = zeros(numFiles,1);
totnumReps = zeros(numFiles,1);
for ii=1:numFiles
    load(fileList(ii).name);
    if exist('MapParams','var') == 1
        totResponse{ii} = MapParams.Response;
        totnumChans(ii) = size(MapParams.significantStimuli,1);
        totnumStimuli(ii) = size(MapParams.significantStimuli,2);
        totnumReps(ii) = size(MapParams.Response,3);

        for jj=1:totnumChans(ii)
            if isnan(MapParams.centerMass.x(jj)) == 0
                xc = MapParams.centerMass.x(jj);
                yc = MapParams.centerMass.y(jj);
                for kk=1:totnumStimuli(ii)
                    xpos = MapParams.centerVals(kk,1);
                    ypos = MapParams.centerVals(kk,2);
                    dist = ceil(sqrt((xc-xpos).^2+(yc-ypos).^2));
                    for ll=1:totnumReps(ii)
                        [minVal,minInd] = min(MapParams.Response(jj,kk,ll,fullWin));
                        if minInd > minWin(1) && minInd < minWin(2) && minVal < -150
                            Pr(dist,1) = Pr(dist,1)+1;
                        end
                        Pr(dist,2) = Pr(dist,2)+1;
                    end
                end
            end
        end
        clear MapParams;
    end
end

distVals = 1:3000;

distVals = distVals(Pr(:,2)>0)';
Pr = Pr(Pr(:,2)>0,:);

distVals = distVals(Pr(:,1)>15);
Pr = Pr(Pr(:,1)>15,:);

finalPr = Pr(:,1)./Pr(:,2);
figure();scatter(distVals,finalPr);

figure();scatter(distVals,log(finalPr));
% Y = A*exp(B*X) ... Pr(VEP-event | distance to center of mass)
% ln(Y) = ln(A)+B*X
%   or Y = A*exp(B*X)+C

N = length(finalPr);

Design = zeros(N,1);
Y = zeros(N,1);
for ii=1:N
    Design(ii) = distVals(ii);
    Y(ii) = -log(finalPr(ii));
end

[b,dev,stats] = glmfit(Design,Y,'normal');


figure();plot(distVals,-log(b(1))*exp(-b(2)*distVals));

% fit non-linear curve
x0 = [1,-0.002,0.2];

myFun = @(x,data) (x(1)*exp(x(2)*data)+x(3));

[x,resnorm] = lsqcurvefit(myFun,x0,distVals,finalPr);

figure();scatter(distVals,finalPr);hold on
plot(1:2000,x(1)*exp(x(2).*(1:2000))+x(3),'r','LineWidth',2);
legend('Data',sprintf('Fit: Pr = %3.3f*exp(%3.3f*x)+%3.3f',x(1),x(2),x(3)));
xlabel('Distance to Retinotopic Center of Mass (pixels)');
ylabel('Probability of VEP-like Event');
title('Likelihood: Pr(VEP-like Event | Distance to Retinotopic Center of Mass)');

x0 = [0.5 200 0.2];

myFun = @(x,data) (x(1)*exp(-(data.*data)./(x(2).*x(2)))+x(3));

[x,resnorm] = lsqcurvefit(myFun,x0,distVals,finalPr);
figure();scatter(distVals,finalPr);hold on;
plot(1:2000,myFun(x,1:2000),'r','LineWidth',2);
legend('Data',sprintf('Fit: Pr = %3.3f*exp(-x^2/(%3.3f)^2)+%3.3f',x(1),x(2),x(3)));
xlabel('Distance to Retinotopic Center of Mass (pixels)');
ylabel('Probability of VEP-like Event');
title('Likelihood: Pr(VEP-like Event | Distance to Retinotopic Center of Mass)');
