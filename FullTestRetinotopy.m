% FullTestRetinotopy.m

% initialize variables
flashes = [250,500,1000,2000];

Center = [500,500];
STDS = [50,150,250,500,1000];
PrHitNoise = 0.2;

PrSeparation = [0.01,0.1,0.2,0.3,0.4,0.5];

numParameters = 4;

hyperParameterFun = @(b,distToCenterMass) (b(1)*exp(-(distToCenterMass.^2)./(2*b(2)*b(2)))+b(3));
DistFun = @(stimCenter,centerVals) (ceil(sqrt((stimCenter(1)-centerVals(:,1)).^2+(stimCenter(2)-centerVals(:,2)).^2))+1);
xaxis = 1:10:2560;
yaxis = 1:10:1440;

lenx = length(xaxis);
leny = length(yaxis);

centerVals = zeros(lenx*leny,2);

count = 1;
for ii=1:lenx
    for jj=1:leny
        centerVals(count,1) = xaxis(ii);
        centerVals(count,2) = yaxis(jj);
        count = count+1;
    end
end
allPossDists = DistFun(Center,centerVals);
maxDist = ceil(sqrt((xaxis(end)-xaxis(1)).^2+(yaxis(end)-yaxis(1)).^2));

stimPrior = ones(lenx*leny,1)./(lenx*leny);

avRelativeError = zeros(length(STDS),length(PrSeparation),length(flashes));
% start looping through different numbers of flashes and different
%   standard deviations and separation parameter values
parfor xx=1:length(STDS)
    display(sprintf('STD = %d',STDS(xx)));
    for yy=1:length(PrSeparation)
        display(sprintf('Separation = %3.2f',PrSeparation(yy)));
        testDist = hyperParameterFun([PrSeparation(yy),STDS(xx),PrHitNoise],(0:maxDist)');
        HitMissProbs = testDist(allPossDists);
        trueParameters = [PrSeparation(yy),Center(1),Center(2),STDS(xx)];
        for zz=1:length(flashes)
            display(sprintf('Number of flashes = %d',flashes(zz)));
            avError = zeros(10,1);
            for avRepeat = 1:5
                % simulate data
                prevIndex = 1;
                Response = zeros(1,flashes(zz),2);
                for ii=1:flashes(zz)
                    dists = DistFun(centerVals(prevIndex,:),centerVals);
                    
                    newPrior = stimPrior;
                    newPrior(dists<100) = 0;
                    newPrior = newPrior./sum(newPrior);
                    unifRand = rand;
                    CDF = cumsum(newPrior(:,1));
                    CDF = CDF-min(CDF);
                    temp = unifRand-CDF;
                    temp(temp<0) = 0;
                    [~,index] = min(temp);
                    prevIndex = index;
                    
                    value = rand;
                    if value <= HitMissProbs(index)
                        hit_or_miss = 1;
                    else
                        hit_or_miss = 0;
                    end
                    Response(1,ii,:) = [index,hit_or_miss];
                end
                
                % calculate maximum likelihood estimate
                [finalParameters] = FitLFPretinoModel(Response,xaxis,yaxis,PrHitNoise,centerVals);
                
                % calculate percent error on the estimate
                error = zeros(numParameters,1);
                for jj=1:numParameters
                    error(jj) = abs(finalParameters(jj)-trueParameters(jj))./trueParameters(jj);
                end
                avError(avRepeat) = mean(error);
            end
            avRelativeError(xx,yy,zz) = mean(avError);
        end
    end
end

save('RetinoModelError.mat','avRelativeError','STDS','PrSeparation','flashes');


versusFlashes = mean(mean(avRelativeError,1),1);
versusSTD = mean(mean(avRelativeError,3),2);
versusSep = mean(mean(avRelativeError,3),1);
figure();subplot(3,1,1);plot(flashes,versusFlashes,'LineWidth',2);
title('Relative Error versus Number of Flashes');
subplot(3,1,2);plot(STDS,versusSTD,'LineWidth',2);
title('Relative Error versus Standard Deviation');
subplot(3,1,3);plot(PrSeparation,versusSep,'LineWidth',2);
title('Relative Error versus Separation Probability Parameter');

acrossSTDS = mean(avRelativeError,1);
figure();imagesc(PrSeparation,flashes,acrossSTDS);
title('Relative Error with STD Averaged Out');

acrossSep = mean(avRelativeError,2);
figure();imagesc(STDS,flashes,acrossSep);
title('Relative Error with Separation Averaged Out');

acrossFlashes = mean(avRelativeError,3);
figure();imagesc(STDS,PrSeparation,acrossFlashes);
title('Relative Error with Flashes Averaged Out');