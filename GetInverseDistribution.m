% get inverse distribution for ProbData to make prior

load('ProbData_RecordingRoomScreen.mat');

InvProbData = zeros(size(ProbData));
numPixels = length(ProbData);
for ii=1:numPixels
    mypdf = ProbData(ii,:);
    newpdf = zeros(length(mypdf),1);
    indeces = find(mypdf~=0);
    mypdf = mypdf(indeces);
    uniformcdf = 1/length(mypdf):1/length(mypdf):1;
%     uniformcdf = uniformcdf(indeces);
    possDists = indeces;
    mycdf = cumsum(mypdf);
    newcdf = (uniformcdf-mycdf)+uniformcdf;
    lastVal = 0;
    for jj=1:length(mycdf)
%         temp = uniformcdf(jj)-mycdf(jj);
%         if temp < 0
%             newcdf(jj) = uniformcdf(jj)+temp;
%         elseif temp == 0
%             newcdf(jj) = uniformcdf(jj);
%         elseif temp > 0
%             newcdf(jj) = uniformcdf(jj)+temp;
%         end
        temp = newcdf(jj)-lastVal;
        newpdf(indeces(jj)) = newcdf(jj)-lastVal;
        lastVal = newcdf(jj);
    end
    temp = newpdf(indeces);
    temp = temp+2.*abs(min(temp));
    temp = temp./sum(temp);
    InvProbData(ii,indeces) = temp';
end

save('InvProbData_RecordingRoomScreen.mat','InvProbData');