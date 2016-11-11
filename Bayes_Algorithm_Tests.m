% Bayes_Algorithm_Tests.m
load('ProbData_RecordingRoomScreen.mat');

b = [0.41 251 .1474];
maxProbHit = b(1)+b(3);

% maybe what I really have is 1000 parameters, which each need to be
% estimated, the probability of a hit at each point in space given the data
% of hits and misses
hyperParameterFun = @(b,distToCenterMass) (b(1)*exp(-(distToCenterMass.^2)./(2*b(2)*b(2)))+b(3));

%likelihoodFun = @(distToCenterMass,SigmaSquare) (repmat((1/sqrt(2*pi)).*SigmaSquare.^(1/2),length(distToCenterMass),1).*exp(-(distToCenterMass.^2)*(SigmaSquare./2)));

DistFun = @(stimCenter,centerVals) (ceil(sqrt((stimCenter(1)-centerVals(:,1)).^2+(stimCenter(2)-centerVals(:,2)).^2))+1);
xaxis = 1:10:2560;
yaxis = 1:10:1440;

lenx = length(xaxis);
leny = length(yaxis);

likelyhit = zeros(lenx,leny);
likelymiss = zeros(lenx,leny);

stimVals = ones(lenx,leny);

xcen = (max(xaxis)-1)/2;
ycen = (max(yaxis)-1)/2;

% the problem with this algorithm is that it creates a map no matter what,
%  it appears to work under conditions when the data is real, but it also 
%  gives what look to be false positives
total = 1;
reps = 5000;
centers = zeros(total,2);
width = zeros(total,2);
alpha = 0.05;

Hitmax = min(maxProbHit,b(3)+b(1));
b(1) = Hitmax-b(3);
maxDist = ceil(sqrt((xaxis(end)-xaxis(1)).^2+(yaxis(end)-yaxis(1)).^2));
pBernoulli = hyperParameterFun(b,(0:maxDist)');

PrHit = mean(pBernoulli);
PrMiss = mean(1-pBernoulli);

% Hitmin = min(pBernoulli);
% Hitmax = max(pBernoulli);
% 
% Missmin = min(1-pBernoulli);
% Missmax = max(1-pBernoulli);
% 
% hitedges = Hitmin-0.01:0.01:Hitmax+0.01;
% missedges = Missmin-0.01:0.01:Missmax+0.01;
% ProbHitUpdate = zeros(size(ProbData));
% ProbMissUpdate = zeros(size(ProbData));
% for ii=1:length(ProbData)
%     temp = ProbData(ii,:)';
%     temp(temp>0) = 1;
%     tempHit = temp.*pBernoulli;
%     tempMiss = temp.*(1-pBernoulli);
%     Hitcounts = histcounts(tempHit,hitedges);
%     Misscounts = histcounts(tempMiss,missedges);
%     Hitcounts = (Hitcounts./sum(Hitcounts));
%     Misscounts = (Misscounts./sum(Misscounts));
%     Hitdiscrete = discretize(tempHit,hitedges);
%     Missdiscrete = discretize(tempMiss,missedges);
%     Hitdiscrete(isnan(Hitdiscrete)) = 1;
%     Missdiscrete(isnan(Missdiscrete)) = 1;
%     newHitProbs = Hitcounts(Hitdiscrete);
%     newMissProbs = Misscounts(Missdiscrete);
% 
%     ProbHitUpdate(ii,:) = newHitProbs;
%     ProbMissUpdate(ii,:) = newMissProbs;
% end


% conditionalProbHit = zeros(size(ProbData));
% conditionalProbMiss = zeros(size(ProbData));
% for ii=1:length(ProbData)
%     tempHit = ProbData(ii,:)'.*pBernoulli;
%     tempMiss = ProbData(ii,:)'.*(1-pBernoulli);
%     figure();plot(tempHit)
%     figure();plot(tempMiss)
%     temp = ProbData(ii,:)';
%     totalProb = sum(tempHit)+sum(tempMiss);
%     conditionalProbHit(ii,:) = pBernoulli./temp; % ./totalProb (which should be one already)
%     conditionalProbMiss(ii,:) = (1-pBernoulli)./temp;
% %     conditionalProbHit(ii,:) = conditionalProbHit(ii,:)./max(conditionalProbMiss(ii,:));
% %     conditionalProbMiss(ii,:) = conditionalProbMiss(ii,:)./max(conditionalProbMiss(ii,:));
% end

centerVals = zeros(lenx*leny,2);
Prior = ones(lenx*leny,1);

count = 1;
for ii=1:lenx
    for jj=1:leny
        centerVals(count,1) = xaxis(ii);
        centerVals(count,2) = yaxis(jj);
%         if ismember(xaxis(ii),1:100) == 0 && ismember(yaxis(jj),1:100) == 0 && ismember(yaxis(jj),yaxis(end)-100:yaxis(end))==0 && ismember(xaxis(ii),xaxis(end)-100:xaxis(end)) == 0
%             Prior(count) = 1;
%         end
        count = count+1;
    end
end

center = [500,600];
allPossDists = DistFun(center,centerVals);
[~,index] = min(allPossDists);
%newP = b(3)*ones(length(pBernoulli),1);
testDist = hyperParameterFun([0.4,250,0.15],(0:maxDist)');
HitMissProbs = testDist(allPossDists);

Prior = Prior./sum(Prior);
stimPrior = ones(lenx*leny,1)./(lenx*leny);
prevIndex = 1;
for zz=1:total
    Response = zeros(reps,3);
    ProbDist = zeros(size(ProbData));
    ProbHitDist = zeros(size(ProbData));
    for ii=1:reps  
        dists = DistFun(centerVals(prevIndex,:),centerVals);
        
        newPrior = Prior;
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
        Response(ii,:) = [centerVals(index,:),hit_or_miss];
        allPossDists = DistFun(Response(ii,1:2),centerVals);
        Inds = (allPossDists-1)*length(allPossDists)+(1:length(allPossDists))';
        ProbDist(Inds) = ProbDist(Inds)+1;
        if hit_or_miss == 1
            ProbHitDist(Inds) = ProbHitDist(Inds)+1;
        end
    end
    for ww=1
        ProbDist = ProbDist./reps;
        ProbHitDist = ProbHitDist./reps;
        
        ProbDist(ProbDist==0) = 1;
        conditionalProb = ProbHitDist./ProbDist;
        
        loglikely = zeros(lenx*leny,1);
        for ii=1:reps
            allPossDists = DistFun(Response(ii,1:2),centerVals);
            Inds = (allPossDists-1)*length(allPossDists)+(1:length(allPossDists))';
            loglikely = loglikely+log(((pBernoulli(allPossDists)).^Response(ii,3)).*((1-pBernoulli(allPossDists)).^(1-Response(ii,3))));
        end

        [B] = PoissonRetinoMap(Response,[max(xaxis),max(yaxis)]);
        
        loglikely = reshape(loglikely,[leny,lenx]);
        figure();subplot(3,1,1);imagesc(xaxis,yaxis,loglikely);set(gca,'YDir','normal');
        
        loglikely2 = zeros(lenx*leny,1);
        for ii=1:100
            inds = randperm(reps);
            temp = Response(inds,3);
            Response(:,3) = temp;
            for jj=1:reps
                allPossDists = DistFun(Response(jj,1:2),centerVals);
                Inds = (allPossDists-1)*length(allPossDists)+(1:length(allPossDists))';
                loglikely2 = loglikely2+log(((pBernoulli(allPossDists)).^Response(jj,3)).*((1-pBernoulli(allPossDists)).^(1-Response(jj,3))));
            end
        end
        
        loglikely2 = loglikely2./100;
        

        loglikely2 = reshape(loglikely2,[leny,lenx]);
        subplot(3,1,2);imagesc(xaxis,yaxis,loglikely2);set(gca,'YDir','normal');
        
        subplot(3,1,3);imagesc(xaxis,yaxis,loglikely-loglikely2);set(gca,'YDir','normal');
        
        unifProb = 1./(lenx*leny);
        Posterior = exp(loglikely-loglikely2);
        Posterior = Posterior./sum(sum(Posterior));
        Posterior = Posterior./unifProb;
        figure();imagesc(xaxis,yaxis,Posterior);set(gca,'YDir','normal');
        caxis([0,50]);
            
%         loglikely = zeros(lenx*leny,1);
%         likely = ones(lenx*leny,1);
%         totalProbData = zeros(lenx*leny,1);
% 
%         for ll=1:reps
%             allPossDists = DistFun(Response(ll,1:2),centerVals);
%             Inds = (allPossDists-1)*length(allPossDists)+(1:length(allPossDists))';
%             loglikely = loglikely+log(((pBernoulli(allPossDists)).^Response(ll,3)).*(((1-pBernoulli(allPossDists))).^(1-Response(ll,3))));
%             loglikely = loglikely+log(((pBernoulli(allPossDists)).^Response(ll,3)).*(((1-pBernoulli(allPossDists))).^(1-Response(ll,3))));
% %             loglikely = loglikely+(((pBernoulli(allPossDists)).^Response(ll,3)).*((1-pBernoulli(allPossDists)).^(1-Response(ll,3))));
% %             totalProbData = totalProbData+((pBernoulli(allPossDists)).^Response(ll,3)).*(((1-pBernoulli(allPossDists))).^(1-Response(ll,3)));
% %             totalProbData = totalProbData+((pBernoulli(allPossDists).*ProbData(Inds)).^Response(ll,3)).*(((1-pBernoulli(allPossDists)).*ProbData(Inds)).^(1-Response(ll,3)));
%             totalProbData = totalProbData+log(((ProbHitUpdate(Inds)).^Response(ll,3)).*((ProbMissUpdate(Inds)).^(1-Response(ll,3))));
% %             loglikely = loglikely+log((((pBernoulli(allPossDists)./(PrHit.*ProbHitUpdate(Inds))).^Response(ll,3)).*(((1-pBernoulli(allPossDists))./(PrMiss.*ProbMissUpdate(Inds))).^(1-Response(ll,3)))));
% %             likely = likely.*(((pBernoulli(allPossDists)./(PrHit)).^Response(ll,3)).*(((1-pBernoulli(allPossDists))./PrMiss).^(1-Response(ll,3))));
% %             loglikely = loglikely+log((((pBernoulli(allPossDists)./(PrHit)).^Response(ll,3)).*(((1-pBernoulli(allPossDists))./PrMiss).^(1-Response(ll,3)))));
% %             loglikely = loglikely+log((((conditionalProbHit(Inds)).^Response(ll,3)).*((conditionalProbMiss(Inds)).^(1-Response(ll,3)))));
% %             likely = likely.*(((conditionalProbHit(Inds)).^Response(ll,3)).*((conditionalProbMiss(Inds)).^(1-Response(ll,3))));
% %             if mod(ll,reps/20) == 0
% %                 z = reshape(likely./totalProbData,[leny,lenx])';
% %                 figure();imagesc(xaxis,yaxis,z');set(gca,'YDir','normal');colorbar;
% %                 title(sprintf('Hit or Miss: %d',Response(ll,3)));
% %             end
%         end
%         Posterior = loglikely-totalProbData;
%         unifProb = 1./(lenx*leny);
%         
%         loglikely = reshape(loglikely,[leny,lenx])';
%         totalProbData = reshape(totalProbData,[leny,lenx])';
%         figure();imagesc(xaxis,yaxis,loglikely');set(gca,'YDir','normal');colorbar;
%         figure();imagesc(xaxis,yaxis,totalProbData');set(gca,'YDir','normal');colorbar;
%         figure();imagesc(xaxis,yaxis,loglikely'-totalProbData');set(gca,'YDir','normal');colorbar;
%         Posterior = exp(Posterior);
%         Posterior = Posterior./sum(Posterior);
%         
%         Posterior = Posterior./unifProb;
%         Posterior = reshape(Posterior,[leny,lenx])';
%         figure();imagesc(xaxis,yaxis,Posterior');set(gca,'YDir','normal');colorbar;
%         caxis([0 100]);
%         
        
        nonZeroInds = find(Response(:,3)==1);
        zeroInds = find(Response(:,3)==0);
        figure();scatter(Response(nonZeroInds,1),Response(nonZeroInds,2),'ob');
        axis([0 max(xaxis) 0 max(yaxis)]);
        hold on;scatter(Response(zeroInds,1),Response(zeroInds,2),'xr');
        
        xdist = sum(Posterior,2);
        ydist = sum(Posterior,1);
        
        halfMaxx = max(xdist)/2;halfMaxy = max(ydist)/2;
        [~,indx] = min(abs(xdist-halfMaxx));
        val1 = xaxis(indx);xdist(indx) = 0;
        [~,indx] = min(abs(xdist-halfMaxx));xdist(indx) = 0;
        val2 = xaxis(indx);fwhm = abs(val2-val1);stan = fwhm/2.355;
        width(zz,1) = atan(((0.236*4*stan)/10)/25)*180/pi;
        
        [~,indy] = min(abs(ydist-halfMaxy));
        val1 = yaxis(indy);ydist(indy) = 0;
        [~,indy] = min(abs(ydist-halfMaxy));ydist(indy) = 0;
        val2 = yaxis(indy);fwhm = abs(val2-val1);stan = fwhm/2.355;
        width(zz,2) = atan(((0.236*4*stan)/10)/25)*180/pi;
        
        [~,centers(zz,1)] = max(sum(Posterior,2));
        [~,centers(zz,2)] = max(sum(Posterior,1));
    end
end
%figure();scatter(xaxis(centers(:,1)),yaxis(centers(:,2)));

