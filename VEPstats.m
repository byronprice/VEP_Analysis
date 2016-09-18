% VEPstats

fileList = dir('RetinoMap*_*.mat');
numFiles = length(fileList);

vepminLats = [];
vepminVals = [];
novepminLats = [];
novepminVals = [];

vepmaxLats = [];
vepmaxVals = [];
novepmaxLats = [];
novepmaxVals = [];
VEP = [];
noVEP = [];

fullWin = 1:200;
minWin = 1:250; % milliseconds
maxWin = 1:250;
totResponse = cell(numFiles,1);
totnumChans = zeros(numFiles,1);
totnumStimuli = zeros(numFiles,1);
totnumReps = zeros(numFiles,1);
totsignifs = cell(numFiles,1);
totminLatency = cell(numFiles,1);
totminVal = cell(numFiles,1);
totmaxLatency = cell(numFiles,1);
totmaxVal = cell(numFiles,1);
for ii=1:numFiles
    load(fileList(ii).name);
    if exist('MapParams','var') == 1
        totResponse{ii} = MapParams.Response;
        totsignifs{ii} = MapParams.significantStimuli > 0;
        totnumChans(ii) = size(MapParams.significantStimuli,1);
        totnumStimuli(ii) = size(MapParams.significantStimuli,2);
        totnumReps(ii) = size(MapParams.Response,3);
        
        totminLatency{ii} = zeros(totnumChans(ii),totnumStimuli(ii),totnumReps(ii));
        totminVal{ii} = zeros(totnumChans(ii),totnumStimuli(ii),totnumReps(ii));
        totmaxLatency{ii} = zeros(totnumChans(ii),totnumStimuli(ii),totnumReps(ii));
        totmaxVal{ii} = zeros(totnumChans(ii),totnumStimuli(ii),totnumReps(ii));
        for jj=1:totnumChans(ii)
            for kk=1:totnumStimuli(ii)
%                 len = size(squeeze(MapParams.Response(jj,kk,:,:)),2);
%                 correlation = corrcoef(squeeze(MapParams.Response(jj,kk,:,:)));
%                 d = -20:20;
%                 B = spdiags(zeros(len,length(d)),d,correlation);
%                 correlation = full(B);
%                 [minvals,ind] = min(correlation(:));
%                 [x,y] = ind2sub(size(correlation),ind);
%                 minlats = abs(x-y);
%                 [maxvals,ind] = max(correlation(:));
%                 [x,y] = ind2sub(size(correlation),ind);
%                   maxlats = abs(x-y)
                [maxvals,maxlats] = max(squeeze(MapParams.Response(jj,kk,:,maxWin)),[],2);
                [minvals,minlats]  = min(squeeze(MapParams.Response(jj,kk,:,minWin)),[],2);
                [maxvals,maxlats] = max(diff(squeeze(MapParams.Response(jj,kk,:,maxWin)),1,2),[],2);
                [minvals,minlats] = min(diff(squeeze(MapParams.Response(jj,kk,:,minWin)),1,2),[],2);
                totminLatency{ii}(jj,kk,:) = minlats+minWin(1)-1;
                totmaxLatency{ii}(jj,kk,:) = maxlats+maxWin(1)-1;
                totmaxVal{ii}(jj,kk,:) = maxvals;
%                 vals = mean(squeeze(MapParams.Response(jj,kk,:,25:50)),2)-min(squeeze(MapParams.Response(jj,kk,:,minWin)),[],2);
                minvals = -minvals;
%                 totminVal{ii}(jj,kk,:) = minvals;
                if totsignifs{ii}(jj,kk) > 0
                    vepminLats = [vepminLats;minlats];
                    vepminVals = [vepminVals;minvals];
                    vepmaxLats = [vepmaxLats;maxlats];
                    vepmaxVals = [vepmaxVals;maxvals];
                    
                    VEP = [VEP;squeeze(MapParams.Response(jj,kk,:,fullWin))];
                else
                    novepminLats = [novepminLats;minlats];
                    novepminVals = [novepminVals;minvals];
                    novepmaxLats = [novepmaxLats;maxlats];
                    novepmaxVals = [novepmaxVals;maxvals];
                    noVEP = [noVEP;squeeze(MapParams.Response(jj,kk,:,fullWin))]; 
                end
            end
        end
        clear MapParams;
    end
end

edges = cell(2,1);
edges{1} = linspace(minWin(1),minWin(end),20); % 20 200 or 1 150
edges{2} = linspace(-100,600,20); %-100 600
figure();hist3([vepminLats,vepminVals],'Edges',edges);
figure();hist3([novepminLats,novepminVals],'Edges',edges);

edges{1} = linspace(maxWin(1),maxWin(end),20); % 100 250
edges{2} = linspace(-100,600,20);
figure();hist3([vepmaxLats,vepmaxVals],'Edges',edges);
figure();hist3([novepmaxLats,novepmaxVals],'Edges',edges);

edges{1} = linspace(minWin(1),minWin(end),20); % 1 150
edges{2} = linspace(maxWin(1),minWin(end),20); % 100 250
figure();hist3([vepminLats,vepmaxLats],'Edges',edges);
figure();hist3([novepminLats,novepmaxLats],'Edges',edges);

inds = vepminLats > 20;
vepminLats = vepminLats(inds);
vepminVals = vepminVals(inds);

inds = novepminLats > 20;
novepminLats = novepminLats(inds);
novepminVals = novepminVals(inds);

inds = vepminLats > 20;
vepmaxLats = vepmaxLats(inds);
vepmaxVals = vepmaxVals(inds);

inds = novepminLats > 20;
novepmaxLats = novepmaxLats(inds);
novepmaxVals = novepmaxVals(inds);

minLats = [vepminLats;novepminLats];
minVals = [vepminVals;novepminVals];
maxLats = [vepmaxLats;novepmaxLats];
maxVals = [vepmaxVals;novepmaxVals];

color = [ones(length(vepminLats),1);zeros(length(novepminLats),1)];
figure();scatter3(minLats,maxLats,minVals,[],color);
% figure();imagesc(cov(VEP));
% figure();imagesc(cov(noVEP));
% figure();imagesc(cov(VEP'));
% figure();imagesc(cov(noVEP'));
% figure();scatter(vepminLats,vepminVals);

%color = [ones(length(vepminLats),1);zeros(length(novepminLats),1)];
% [coeff,score] = pca([vepminLats,vepminVals,vepmaxLats,vepmaxVals]);

% edges{1} = linspace(-2000,2000,20);
% edges{2} = linspace(-2000,2000,20);
% figure();scatter3(score(:,1),score(:,2),score(:,3),[]);
% figure();scatter(score(:,1),score(:,2),[]);
% coeff
% 
% figure();scatter3(vepminLats,vepminVals,vepmaxLats,[]);
% coeff

% [coeff,score] = pca(VEP);
% figure();scatter(score(:,1),score(:,2));

PrMinVEPlats = (sum(vepminLats < 110)-sum(vepminLats < 40))./length(vepminLats);
PrMinVEPvals = sum(vepminVals > 50)./length(vepminVals);

inds = vepminLats > 40;
newinds = vepminLats(inds) < 110;

Prconditional = sum(vepminVals(newinds)>50)./length(newinds);

combined = vepminLats(vepminVals > 50);
Prboth = (sum(combined < 110)-sum(combined < 40))./length(vepminLats);



PrMinNoVEPlats = (sum(novepminLats < 110)-sum(novepminLats < 40))./length(novepminLats);
PrMinNoVEPvals = sum(novepminVals > 50)./length(novepminVals);

inds = novepminLats > 40;
newinds = novepminLats(inds) < 110;

Prconditional = sum(novepminVals(newinds)>50)./length(newinds);

combined = novepminLats(novepminVals > 50);
Prboth = (sum(combined < 110)-sum(combined < 40))./length(novepminLats);

window = 60:120;

thresh = mean(2*1.4826*mad(noVEP,1,2));
numCrossings = 0;
for ii=1:size(VEP,1)
    temp = VEP(ii,window)+thresh;
    maxs = max(temp,0);maxs = min(maxs,1);
    mins = min(temp,0);mins = max(mins,-1);
    tot = maxs+mins;
    newTemp = sum(diff(tot) <=-2);
    if newTemp >= 1
        numCrossings = numCrossings+1;
    end
end
Pr = numCrossings/size(VEP,1)

numCrossings = 0;
for ii=1:size(noVEP,1)
    temp = noVEP(ii,window)+thresh;
    maxs = max(temp,0);maxs = min(maxs,1);
    mins = min(temp,0);mins = max(mins,-1);
    tot = maxs+mins;
    newTemp = sum(diff(tot) <=-2);
    if newTemp >= 1
        numCrossings = numCrossings+1;
    end
end
Pr = numCrossings/size(noVEP,1)

