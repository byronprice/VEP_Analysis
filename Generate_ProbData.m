% Generate_ProbData.m

% by either a convolution algorithm or Monte Carlo draws from the
%  distributions, create the matrix that represents the probability of the
%  data for my Bayesian LFP summation field mapping experiments
xaxis = 1:10:2560;
yaxis = 1:10:1440;

xlen = length(xaxis);
ylen = length(yaxis);

screen = zeros(2560,1440);
screen(xaxis,yaxis) = 1;
numPositions = xlen*ylen;
maxdist = ceil(sqrt((xaxis(end)-1).^2+(yaxis(end)-1).^2));

% DistFun = @(center,allPoss) (ceil(sqrt((center-allPoss(:,1)).^2+(center-allPoss(:,2)).^2)));
% ProbData = zeros(numPositions,maxdist+1);
% ProbData(:,1) = 1;
% matSize = 3:2:5847;
% for ii=2:maxdist
%     dist = ii-1;
%     kernel = zeros(matSize(ii-1));
%     center = ceil(matSize(ii-1)/2);
%     
%     allPoss = zeros(matSize(ii-1)*matSize(ii-1),2);
%     count = 1;
%     for kk=1:matSize(ii-1)
%         for jj=1:matSize(ii-1)
%             allPoss(count,:) = [jj,kk];
%             count = count+1;
%         end
%     end
%     dists = DistFun(center,allPoss);
%     kernel(dists==dist) = 1;
%     
%     C = conv2(screen,kernel,'same');
%     C = C(xaxis,yaxis);
%     ProbData(:,ii) = reshape(C',[numPositions,1]);
% end
% 
% for ii=1:numPositions
%     ProbData(ii,:) = ProbData(ii,:)./sum(ProbData(ii,:));
% end
% 
% figure();plot(ProbData(1,:));
% 
% figure();plot(ProbData(2,:));
% figure();plot(ProbData(10,:));
% 
% figure();plot(ProbData(1000,:));
% 
% figure();plot(ProbData(5000,:));

centerVals = zeros(xlen*ylen,2);

count = 1;
for ii=1:xlen
    for jj=1:ylen
        centerVals(count,1) = xaxis(ii);
        centerVals(count,2) = yaxis(jj);
        count = count+1;
    end
end

stimSelection = ones(xlen*ylen,1)./(xlen*ylen);
prevIndex = 1;

% y = randperm(100); % if y is a vector of indeces, this is one way to add
%  elements to the result
% result = zeros(100,100);
% for ii=1:100
% result(ii,y(ii)) = result(ii,y(ii))+1;
% end

numPositions = length(centerVals);
DistConversion = @(Distance) ((Distance(:)-1).*(numPositions)+(1:numPositions)');

DistFun = @(center,allPoss) (ceil(sqrt((center(1)-allPoss(:,1)).^2+(center(2)-allPoss(:,2)).^2))+1);
ProbData = zeros(numPositions,maxdist+1);
reps = 1000;
for ii=1:reps
    dists = DistFun(centerVals(prevIndex,:),centerVals);
    tempInds = DistConversion(dists);
    ProbData(tempInds) = ProbData(tempInds)+1;
    newPrior = stimSelection;
    newPrior(dists<150) = 0;
    newPrior = newPrior./sum(newPrior);
    unifRand = rand;
    CDF = cumsum(newPrior(:,1));
    CDF = CDF-min(CDF);
    temp = unifRand-CDF;
    temp(temp<0) = 0;
    [~,index] = min(temp);
    prevIndex = index;
end
ProbData = ProbData./reps;

for ii=1:100:1000
figure();
plot(ProbData(ii,:));
end