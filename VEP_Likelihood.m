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

lowMin=[60];
highMin =[120];
thresh = [-150];
maxParam = zeros(4,4);
minParam = zeros(4,4);
Residuals = zeros(4,4);

load('VEP_Prototype');

total = [];
for lowii = 1
    for highii = 1
        for threshii=1
            Pr_Hit = zeros(3000,2);
            
            Pr_Size = zeros(3000,2);
            
            sizes = [];
            totSignifs = cell(numFiles,1);
            totResponse = cell(numFiles,1);
            totnumChans = zeros(numFiles,1);
            totnumStimuli = zeros(numFiles,1);
            totnumReps = zeros(numFiles,1);
            for ii=48
                load(fileList(ii).name);
                if exist('MapParams','var') == 1
                    totResponse{ii} = MapParams.Response;
                    totnumChans(ii) = size(MapParams.significantStimuli,1);
                    totnumStimuli(ii) = size(MapParams.significantStimuli,2);
                    totnumReps(ii) = size(MapParams.Response,3);
                    totSignifs{ii} = MapParams.significantStimuli;
                    
                    for jj=1:totnumChans(ii)
                        if isnan(MapParams.centerMass.x(jj)) == 0
                            xc = MapParams.centerMass.x(jj);
                            yc = MapParams.centerMass.y(jj);
                            Y = [];
                            X = [];
                            DISTS = [];
                            for kk=1:totnumStimuli(ii)
                                xpos = MapParams.centerVals(kk,1);
                                ypos = MapParams.centerVals(kk,2);
                                dist = ceil(sqrt((xc-xpos).^2+(yc-ypos).^2));
                                dist = dist-mod(dist,20)+1;
                                
                                for ll=1:totnumReps(ii)
                                    VEP = squeeze(MapParams.Response(jj,kk,ll,fullWin));
                                    baseline = 0;
                                    VEP = VEP-baseline;
                                    
                                    % get scalar multiple for vector
                                    %  division
                                    Y = [Y;dot(VEP,VEP_Prototype)/(norm(VEP_Prototype).^2)];
                                    X = [X;xpos,ypos];
                                    DISTS = [DISTS;dist-1];
                                    
                                    [minVal,minInd] = min(VEP);
                                                                   
%                                     sizes = [sizes;[minVal,dist]];
                                    if minInd > lowMin(lowii) && minInd < highMin(highii) && minVal < thresh(threshii)
                                       Pr_Hit(dist,1) = Pr_Hit(dist,1)+1;
                                       total = [total;[1,dist-1]];
                                    else
                                       total = [total;[0,dist-1]]; 
                                    end
                                    Pr_Hit(dist,2) = Pr_Hit(dist,2)+1;
                                    if minVal < thresh(threshii)
                                        Pr_Size(dist,1) = Pr_Size(dist,1)+1;
                                    end
                                    Pr_Size(dist,2) = Pr_Size(dist,2)+1;
                                end
                            end
                        
                            % Uri's VEP model non-linear regression
                            b0 = [700,1200,20000,50,25000];
                            figure();scatter(DISTS,Y);
                            myFun = @(b,x) exp(-0.5.*((b(5)/(b(3)*b(5)-b(4)*b(4))).*(x(:,1)-b(1)).^2-...
                                (x(:,1)-b(1)).*(x(:,2)-b(2)).*(2*b(4)/(b(3)*b(5)-b(4)*b(4)))+...
                                (b(3)/(b(3)*b(5)-b(4)*b(4))).*(x(:,2)-b(2)).^2));
                            
                            [beta,R,J,covB,MSE] = nlinfit(X,Y,myFun,b0);
                            
                            %minFun = @(b) 
                            %fmincon for constrained minimization
                            
                            ci = nlparci(beta,R,'covar',covB);
                            
                            x = 0:10:2500;
                            y = 0:10:1600;
                            Z = zeros(length(x),length(y));
                            for ww=1:length(x)
                                for xx=1:length(y)
                                    Z(ww,xx) = myFun(beta,[x(ww),y(xx)]);
                                end
                            end
                            figure();imagesc(x,y,Z');set(gca,'YDir','normal');
                        
                        end
                    end
                    clear MapParams;
                end
            end
            
            distVals = 1:3000;
            
            % figure();scatter(distVals,Pr_Size(:,1)./Pr_Size(:,2));
            % title('Probability of VEP < -150 versus Distance to Center of Mass')
            %
            % figure();scatter(sizes(:,2),sizes(:,1));
            % title('Size of VEP versus Distance to Center of Mass')
            
            distVals = distVals(Pr_Hit(:,2)>0)';
            Pr_Hit = Pr_Hit(Pr_Hit(:,2)>0,:);
            
            distVals = distVals(Pr_Hit(:,1)>1);
            Pr_Hit = Pr_Hit(Pr_Hit(:,1)>1,:);
            
            finalPr = Pr_Hit(:,1)./Pr_Hit(:,2);
            % figure();scatter(distVals,finalPr);
            % title('Probability of VEP-like Event versus Distance to Center of Mass')
            %
            % % Y = A*exp(B*X) ... Pr(VEP-event | distance to center of mass)
            % % ln(Y) = ln(A)+B*X
            % %   or Y = A*exp(B*X)+C
            %
            % N = length(finalPr);
            %
            % Design = zeros(N,1);
            % Y = zeros(N,1);
            % for ii=1:N
            %     Design(ii) = distVals(ii);
            %     Y(ii) = -log(finalPr(ii));
            % end
            %
            % [b,dev,stats] = glmfit(Design,Y,'normal');
            %
            %
            % figure();plot(distVals,-log(b(1))*exp(-b(2)*distVals));
            % title(sprintf('Exponential Fit to Probability Data, Pr = %3.2f*exp(%3.2f*X)',-log(b(1)),-b(2)));
            %
            % % fit non-linear curve
            % x0 = [1,-0.002,0.2];
            %
            % myFun = @(x,data) (x(1)*exp(x(2)*data)+x(3));
            %
            % [x,resnorm] = lsqcurvefit(myFun,x0,distVals,finalPr);
            % resnorm
            %
            % figure();scatter(distVals,finalPr);hold on
            % plot(1:2000,x(1)*exp(x(2).*(1:2000))+x(3),'r','LineWidth',2);
            % legend('Data',sprintf('Fit: Pr = %3.3f*exp(%3.3f*x)+%3.3f',x(1),x(2),x(3)));
            % xlabel('Distance to Retinotopic Center of Mass (pixels)');
            % ylabel('Probability of VEP-like Event');
            % title('Likelihood: Pr(VEP-like Event | Distance to Retinotopic Center of Mass)');
            
%             x0 = [0.5 200 0.2];
%             
%             myFun = @(x,data) (x(1)*exp(-(data.*data)./(x(2).*x(2)))+x(3));
%             
%             [x,resnorm] = lsqcurvefit(myFun,x0,distVals,finalPr);
%             resnorm
                
              %x0 = [0.2 0.02 300 0.09];
              x0 = [0.3 250 0.1];
                  
              %myFun = @(x,data) (x(1)./(1+exp(x(2)*(data-x(3))))+x(4));
              myFun = @(x,data) (x(1).*exp(-(data.*data)./(2.*x(2).*x(2)))+x(3));

              [x,resnorm] =  lsqcurvefit(myFun,x0,distVals,finalPr);

            maxParam(lowii,highii,threshii) = x(1);
            minParam(lowii,highii,threshii) = x(3);
            Residuals(lowii,highii,threshii) = resnorm;
            figure();scatter(distVals,finalPr);hold on;
            plot(1:2000,myFun(x,1:2000),'r','LineWidth',2);
            legend('Data',sprintf('Min Latency = %d, Max Latency = %d (msecs)',lowMin(lowii),highMin(highii)));
            xlabel('Distance to Retinotopic Center of Mass (pixels)');
            ylabel('Probability of VEP-like Event');
            title('Likelihood: Pr(VEP-like Event at Specific Distance to Retinotopic Center of Mass)');
        end
     end
end

newTotal = total(total(:,2)<=2000,:);

totalDataPoints = length(newTotal);

allDists = unique(newTotal(:,2));

distLen = length(allDists);

distNums = zeros(distLen,1);

hitAndDistNums = zeros(distLen,1);
missAndDistNums = zeros(distLen,1);

forDistFit = [];
forHitDistFit = [];
forMissDistFit = [];
for ii=1:totalDataPoints
    for jj=1:distLen
        if newTotal(ii,2) == allDists(jj)
            distNums(jj) = distNums(jj)+1;
            forDistFit = [forDistFit,allDists(jj)];
        end
        if newTotal(ii,2) == allDists(jj) && newTotal(ii,1) == 1
            hitAndDistNums(jj) = hitAndDistNums(jj)+1;
            forHitDistFit = [forHitDistFit,allDists(jj)];
        end
        if newTotal(ii,2) == allDists(jj) && newTotal(ii,1) == 0
           missAndDistNums(jj) = missAndDistNums(jj)+1;
        end
    end
end
PrDist = distNums./totalDataPoints;
PrHitAndDist = hitAndDistNums./totalDataPoints;

PrMissAndDist = missAndDistNums./totalDataPoints;


conditionalProb = PrHitAndDist./PrDist;

figure();plot(allDists,PrDist);title('Probability of Being Specific Distance to Center of Mass');

figure();plot(allDists,PrHitAndDist);title('Probability of Hit and Being Specific Distance to Center of Mass');

figure();plot(allDists,PrMissAndDist);title('Probability of Miss and Being Specific Distance to Center of Mass');

% x0 = 500;
%                   
% %myFun = @(x,data) (x(1)./(1+exp(x(2)*(data-x(3))))+x(4));
% myFun = @(x,data) ((x.^(data)*exp(-x))./(factorial(x)));
% 
% [x,resnorm] =  lsqcurvefit(myFun,x0,allDists,PrHitAndDist);

% figure();plot(allDists,myFun(x,allDists));
figure();plot(allDists,myFun(x,allDists).*PrDist);
figure();plot(allDists,(1-myFun(x,allDists)).*PrDist);
figure();plot(allDists,conditionalProb);title('Conditional Probability of Hit | Distance to Center of Mass');
figure();plot(allDists,conditionalProb./sum(conditionalProb));title('Normalized Conditional Probability');

