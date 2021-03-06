% Bayes_Algorithm_Tests.m
%load('InvProbData_RecordingRoomScreen.mat');
load('ProbData_RecordingRoomScreen.mat');

reps = 1000;

trueCenter = [1000,700];
trueSTD = 250;
display(sprintf('True Center Point: %d, %d (screen position in pixels)',trueCenter(1),trueCenter(2)));
display(sprintf('True Standard Deviation: %d pixels',trueSTD));
%likelihoodFun = @(distToCenterMass,SigmaSquare) (repmat((1/sqrt(2*pi)).*SigmaSquare.^(1/2),length(distToCenterMass),1).*exp(-(distToCenterMass.^2)*(SigmaSquare./2)));

DistFun = @(stimCenter,centerVals) (ceil(sqrt((stimCenter(1)-centerVals(:,1)).^2+(stimCenter(2)-centerVals(:,2)).^2))+1);
xaxis = 1:10:2560;
yaxis = 1:10:1440;

lenx = length(xaxis);
leny = length(yaxis);

STDS = 10:50:1000;
lenstd = length(STDS);
h_dist = 10;
h_std = 50;
FullLikelihood = zeros(leny,lenx,lenstd);

bigcount = 1;
for std=STDS
b = [0.4 std .2];
maxProbHit = b(1)+b(3);
PrHitNoise = b(3);
hyperParameterFun = @(b,distToCenterMass) (b(1)*exp(-(distToCenterMass.^2)./(2*b(2)*b(2)))+b(3));


likelyhit = zeros(lenx,leny);
likelymiss = zeros(lenx,leny);

stimVals = ones(lenx,leny);

xcen = (max(xaxis)-1)/2;
ycen = (max(yaxis)-1)/2;

% the problem with this algorithm is that it creates a map no matter what,
%  it appears to work under conditions when the data is real, but it also 
%  gives what look to be false positives
total = 1;

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

center = trueCenter;
allPossDists = DistFun(center,centerVals);
[~,index] = min(allPossDists);
%newP = b(3)*ones(length(pBernoulli),1);
testDist = hyperParameterFun([0.4,trueSTD,PrHitNoise],(0:maxDist)');
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
        
        numParameters = 4;
        
        numRepeats = 10;
        gradientVec = zeros(numParameters,1);
        bigParameterVec = zeros(numRepeats,numParameters);
        maxITER = 1000;
        tolerance = 0.01;
        h = [0.01,1,1,1];
        
        Bounds = [0,1-PrHitNoise;min(xaxis),max(xaxis);min(yaxis),max(yaxis);1,1000];
   
        display('Steepest Ascent ...');
        flashPoints = Response(:,1:2);
        hitMiss = Response(:,3);
        
        % repeat gradient ascent from a number of different starting
        % positions
        for repeats = 1:numRepeats
            parameterVec = zeros(maxITER,numParameters);
            logLikelihood = zeros(maxITER,1);
            %parameterVec(1,:) = squeeze(bigParameterVec(repeats,:));
            parameterVec(1,1) = Bounds(1,1)+(Bounds(1,2)-Bounds(1,1)).*rand;
            parameterVec(1,2) = Bounds(2,1)+(Bounds(2,2)-Bounds(2,1)).*rand;
            parameterVec(1,3) = Bounds(3,1)+(Bounds(3,2)-Bounds(3,1)).*rand;
            parameterVec(1,4) = Bounds(4,1)+(Bounds(4,2)-Bounds(4,1)).*rand;
            
            figure();scatter3(trueCenter(1),trueCenter(2),trueSTD,'filled');axis([0 max(xaxis) 0 max(yaxis) 0 Bounds(4,2)]);hold on;
            xlabel('Horizontal Screen Position');ylabel('Vertical Screen Position');zlabel('Standard Deviation');
            pause(5);
%             figure();plot(b(1),0.75,'*','MarkerSize',20);axis([0 1 0 1]);hold on;
            % for each starting position, do maxITER iterations
            for iter=2:maxITER
                % calcualte likelihood at the current position in parameter
                %  space
                for kk=1:reps
                    distance = DistFun(flashPoints(kk,:),[parameterVec(iter-1,2),parameterVec(iter-1,3)]);
                    p = hyperParameterFun([parameterVec(iter-1,1),parameterVec(iter-1,4),PrHitNoise],distance);
                    logLikelihood(iter-1) = logLikelihood(iter-1)+log((p.^(hitMiss(kk))).*((1-p).^(1-hitMiss(kk))));
                end
                
                if iter > 300
                    check = sum(diff(logLikelihood(iter-250:iter-1)));
                    if check < tolerance
                        break;
                    end
                end
                
                temp = parameterVec(iter-1,:);
                randInd = random('Discrete Uniform',numParameters,1);
                temp(randInd) = Bounds(randInd,1)+(Bounds(randInd,2)-Bounds(randInd,1)).*rand;
%                 for kk=1:numParameters
%                    temp(kk) = Bounds(kk,1)+(Bounds(kk,2)-Bounds(kk,1)).*rand; 
%                 end
                likely = 0;
                for kk=1:reps
                    distance = DistFun(flashPoints(kk,:),[temp(2),temp(3)]);
                    p = hyperParameterFun([temp(1),temp(4),PrHitNoise],distance);
                    likely = likely+log((p.^(hitMiss(kk))).*((1-p).^(1-hitMiss(kk))));
                end
                
                if likely > logLikelihood(iter-1)
                    parameterVec(iter-1,:) = temp';
                    logLikelihood(iter-1) = likely;
                end
                
            
                    % calculate the gradient by calculating the likelihood
                    %  after moving over a small step h along each
                    %  dimension
                    for jj=1:numParameters
                        tempParameterVec = parameterVec(iter-1,:);
                        tempParameterVec(jj) = tempParameterVec(jj)+h(jj);
                        gradLikelihood = 0;
                        for kk=1:reps
                            distance = DistFun(flashPoints(kk,:),[tempParameterVec(2),tempParameterVec(3)]);
                            p = hyperParameterFun([tempParameterVec(1),tempParameterVec(4),PrHitNoise],distance);
                            gradLikelihood = gradLikelihood+log((p.^(hitMiss(kk))).*((1-p).^(1-hitMiss(kk))));
                        end
                        gradientVec(jj) = (gradLikelihood-logLikelihood(iter-1))./h(jj);
                    end
                    
                    % line search to get distance to move along gradient
                    alpha = [0,1e-6,1e-4,1e-2,1e-1,1e0,1e1];
                    lineSearchLikelihoods = zeros(length(alpha),1);
                    lineSearchLikelihoods(1) = logLikelihood(iter-1);
                    
                    for ii=2:length(alpha)
                        tempParameterVec = parameterVec(iter-1,:)'+gradientVec.*alpha(ii);
                        for kk=1:reps
                            distance = DistFun(flashPoints(kk,:),[tempParameterVec(2),tempParameterVec(3)]);
                            p = hyperParameterFun([tempParameterVec(1),tempParameterVec(4),PrHitNoise],distance);
                            lineSearchLikelihoods(ii) = lineSearchLikelihoods(ii)+log((p.^(hitMiss(kk))).*((1-p).^(1-hitMiss(kk))));
                        end
                    end
                    
                    newInds = find(imag(lineSearchLikelihoods)==0);
                    newValues = zeros(length(newInds),1);
                    for ii=1:length(newInds)
                        newValues(ii) = lineSearchLikelihoods(newInds(ii));
                    end
                    [~,ind] = max(newValues);
                    ind = newInds(ind);
                    
                    for jj=1:numParameters
                        parameterVec(iter,jj) = max(Bounds(jj,1),min(parameterVec(iter-1,jj)+alpha(ind)*gradientVec(jj),Bounds(jj,2))); 
                    end
                    title(sprintf('b_1 = %3.2f',parameterVec(iter,1)));
                    scatter3(parameterVec(iter,2),parameterVec(iter,3),parameterVec(iter,4),'filled');pause(0.05);
%                     scatter(parameterVec(iter,1),0.5,'filled');pause(0.05);
            end
            bigParameterVec(repeats,:) = parameterVec(iter-1,:);
        end

        median(bigParameterVec,1)
%         loglikely = zeros(lenx*leny,1);
%         loginv = zeros(lenx*leny,1);
%         logprobs = zeros(lenx*leny,1);
%         for ii=1:reps
%             allPossDists = DistFun(Response(ii,1:2),centerVals);
%             Inds = (allPossDists-1)*length(allPossDists)+(1:length(allPossDists))';
%             loginv = loginv+log(InvProbData(Inds));
%             logprobs = logprobs+log(ProbData(Inds));
%             loglikely = loglikely+log(((pBernoulli(allPossDists)).^Response(ii,3))...
%                 .*(((1-pBernoulli(allPossDists))).^(1-Response(ii,3))));
%         end
%         
%  %       [B] = PoissonRetinoMap(Response,[max(xaxis),max(yaxis)]);
%         loginv = reshape(loginv,[leny,lenx]);
%         logprobs = reshape(logprobs,[leny,lenx]);
%         
%         loglikely = reshape(loglikely,[leny,lenx]);
%         
%         FullLikelihood(:,:,bigcount) = loglikely;
%         bigcount = bigcount+1;
 %       figure();imagesc(xaxis,yaxis,loglikely);set(gca,'YDir','normal');
   %     colorbar;
%         figure();imagesc(xaxis,yaxis,loginv);set(gca,'YDir','normal');colorbar;
%         figure();imagesc(xaxis,yaxis,logprobs);set(gca,'YDir','normal');colorbar;
        
%         loglikely2 = zeros(lenx*leny,1);
%         for ii=1:100
%             inds = randperm(reps);
%             temp = Response(inds,3);
%             Response(:,3) = temp;
%             for jj=1:reps
%                 allPossDists = DistFun(Response(jj,1:2),centerVals);
%                 Inds = (allPossDists-1)*length(allPossDists)+(1:length(allPossDists))';
%                 loglikely2 = loglikely2+log(((pBernoulli(allPossDists).*InvProbData(Inds)).^Response(ii,3))...
%                 .*(((1-pBernoulli(allPossDists)).*InvProbData(Inds)).^(1-Response(ii,3))));
%             end
%         end
%         
%         loglikely2 = loglikely2./100;
%         
% 
%         loglikely2 = reshape(loglikely2,[leny,lenx]);
%         subplot(3,1,2);imagesc(xaxis,yaxis,loglikely2);set(gca,'YDir','normal');
%         colorbar;
%         
%         subplot(3,1,3);imagesc(xaxis,yaxis,loglikely-loglikely2);set(gca,'YDir','normal');
%         colorbar;
%         
%         unifProb = 1./(lenx*leny);
%         Posterior = exp(loglikely-loglikely2);
%         Posterior = Posterior./sum(sum(Posterior));
%         Posterior = Posterior./unifProb;
%         figure();imagesc(xaxis,yaxis,Posterior);set(gca,'YDir','normal');
%         caxis([0,50]);
            
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
        
%         nonZeroInds = find(Response(:,3)==1);
%         zeroInds = find(Response(:,3)==0);
%         figure();scatter(Response(nonZeroInds,1),Response(nonZeroInds,2),'ob');
%         axis([0 max(xaxis) 0 max(yaxis)]);
%         hold on;scatter(Response(zeroInds,1),Response(zeroInds,2),'xr');
        
%         xdist = sum(Posterior,2);
%         ydist = sum(Posterior,1);
%         
%         halfMaxx = max(xdist)/2;halfMaxy = max(ydist)/2;
%         [~,indx] = min(abs(xdist-halfMaxx));
%         val1 = xaxis(indx);xdist(indx) = 0;
%         [~,indx] = min(abs(xdist-halfMaxx));xdist(indx) = 0;
%         val2 = xaxis(indx);fwhm = abs(val2-val1);stan = fwhm/2.355;
%         width(zz,1) = atan(((0.236*4*stan)/10)/25)*180/pi;
%         
%         [~,indy] = min(abs(ydist-halfMaxy));
%         val1 = yaxis(indy);ydist(indy) = 0;
%         [~,indy] = min(abs(ydist-halfMaxy));ydist(indy) = 0;
%         val2 = yaxis(indy);fwhm = abs(val2-val1);stan = fwhm/2.355;
%         width(zz,2) = atan(((0.236*4*stan)/10)/25)*180/pi;
%         
%         [~,centers(zz,1)] = max(sum(Posterior,2));
%         [~,centers(zz,2)] = max(sum(Posterior,1));
    end
end
%figure();scatter(xaxis(centers(:,1)),yaxis(centers(:,2)));
end

[M,I] = max(FullLikelihood(:));
[I,J,K] = ind2sub(size(FullLikelihood),I);
 
display(M);
display(I);display(J);display(K);

reducedLikely = FullLikelihood(I-1:I+1,J-1:J+1,K-1:K+1);
FisherInfo = zeros(3,3);
% [ sigma_xx sigma_xy sigma_xstd        ]
% [ sigma_yx sigma_yy sigma_ystd        ] 
% [ sigma_stdx sigma_stdy sigma_stdstd  ]

FisherInfo(1,1) = (reducedLikely(2,2+1,2)-2*reducedLikely(2,2,2)+...
    reducedLikely(2,2-1,2))./(h_dist^2);
FisherInfo(2,2) = (reducedLikely(2+1,2,2)-2*reducedLikely(2,2,2)+...
    reducedLikely(2-1,2,2))./(h_dist^2);
FisherInfo(3,3) = (reducedLikely(2,2,2+1)-2*reducedLikely(2,2,2)+...
    reducedLikely(2,2,2-1))./(h_std^2);
FisherInfo(1,2) = (reducedLikely(2+1,2+1,2)-reducedLikely(2+1,2-1,2)-...
    reducedLikely(2-1,2+1,2)+reducedLikely(2-1,2-1,2))./(4*h_dist^2);
FisherInfo(2,1) = FisherInfo(1,2);
FisherInfo(1,3) = (reducedLikely(2,2+1,2+1)-reducedLikely(2,2+1,2-1)-...
    reducedLikely(2,2-1,2+1)+reducedLikely(2,2-1,2-1))./(4*h_dist*h_std);
FisherInfo(3,1) = FisherInfo(1,3);
FisherInfo(2,3) = (reducedLikely(2+1,2,2+1)-reducedLikely(2+1,2,2-1)-...
    reducedLikely(2-1,2,2+1)+reducedLikely(2-1,2,2-1))./(4*h_dist*h_std);
FisherInfo(3,2) = FisherInfo(2,3);

covMat = inv(-FisherInfo);

maxX = xaxis(J);maxY = yaxis(I);maxSTD = STDS(K);
errorX = sqrt(covMat(1,1));errorY = sqrt(covMat(2,2));
errorSTD = sqrt(covMat(3,3));
display(sprintf('Estimate for X position: %4.3f +/- %4.3f',maxX,errorX.*1.96));
display(sprintf('Estimate for Y position: %4.3f +/- %4.3f',maxY,errorY.*1.96));
display(sprintf('Estimate for Standard Deviation: %4.3f +/- %4.3f',maxSTD,errorSTD.*1.96));
