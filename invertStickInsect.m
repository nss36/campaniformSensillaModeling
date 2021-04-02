% invertStickInsect
% Nicholas Szczecinski
% Department of Mechanical and Aerospace Engineering
% West Virginia University
% 2 April 2021

clear
close all
clc

load('stickInsectWalkingForces')

load('stick6AParams');
tau6A = stick6AParams.tau;
a6A = stick6AParams.a;
b6A = stick6AParams.b;
c6A = stick6AParams.c;
d6A = stick6AParams.d;

load('stick6BParams');
tau6B = stick6BParams.tau;
a6B = stick6BParams.a;
b6B = stick6BParams.b;
c6B = stick6BParams.c;
d6B = stick6BParams.d;

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ONE-GROUP INVERSION METHOD
% FIGURE 9
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

x01 = 0;
x02 = 0;
k = 0;
lb = -.5;
ub = 3;
p0 = .1;
xTol = 1e-4;
fTol = 1e-10;

for i=5:8    
    hOneGroupFig = figure;
    
    tSamp = cleanedData{i}.t(1:end-1);
    u6A = cleanedData{i}.u(2:end);
    u6B = -cleanedData{i}.u(2:end);
    
    tTW = cleanedData{i}.t(1:end-1);
    yTW = cleanedData{i}.sens6B(2:end);
    yTW(yTW < 10) = 0;
    
    nSamps = length(tSamp);
    uSamp = zeros(nSamps,2);
    
    x0perm = u6B(1) - (yTW(1) - c6B*u6B(1) - d6B)/a6B;
    x01 = x0perm;
    x02 = x0perm;
    
    for j=2:nSamps

        %minimize the force
        x0 = uSamp(j-1,1);
        jmax = j;
        tLoop = tTW(1:jmax);
        yDesTW = yTW(1:jmax);

        g1 = @(x)loadMatchingFunc(x,j,tau6B,a6B,b6B,c6B,d6B,tSamp,uSamp(:,1),yDesTW,x01);

        if g1(x0) == 0
            %do nothing
        else
            x0 = fzero(g1,[lb,ub]);
        end

        [uSamp(j,1),~,~,converged] = cornerSearch(g1,x0,p0,lb,ub,100,xTol,fTol,false);
        if ~converged
            keyboard
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%
        %maximize the force
        
        x0 = uSamp(j-1,2);

        g2 = @(x)loadMatchingFunc(x,j,tau6B,a6B,b6B,c6B,d6B,tSamp,uSamp(:,2),yDesTW,x02);

        if g2(x0) == 0
            %do nothing
        else
            x0 = fzero(g2,[lb,ub]);
        end

        [uSamp(j,2),~,~,converged] = cornerSearch(g2,x0,p0,lb,ub,100,xTol,fTol,true);
        if ~converged
            keyboard
        end

        figure(hOneGroupFig)
        subplot(1,2,1)
        cla
        plot(tLoop,yDesTW,'linewidth',2)
        hold on
        [yLo,xLo] = simulateMinusLowpassPL(uSamp(1:j,1),tau6B,tSamp(1:j),a6B,b6B,c6B,d6B,x0perm);
        plot(tLoop,yLo,':','linewidth',2)
        [yHi,xHi] = simulateMinusLowpassPL(uSamp(1:j,2),tau6B,tSamp(1:j),a6B,b6B,c6B,d6B,x0perm);

        x01 = xLo(end);
        x02 = xHi(end);

        subplot(1,2,2)
        cla
        plot(tLoop,u6B(1:j),'linewidth',2)
        hold on
        plot(tLoop,uSamp(1:j,1),'--','linewidth',2)
        plot(tLoop,uSamp(1:j,2),':','linewidth',2)
        drawnow
        
    end

    uInterped = @(u) pchip(tSamp,u,tTW);

    lineColors = lines(3);

    subplot(1,2,2)
    hold on
    legend('Actual force','Lower bound force','Upper bound force','location','northwest','FontSize',8)
    title('(b) Estimated and actual stimulus force','FontSize',8)
    xlabel('time (s)','FontSize',8)
    ylabel('mN','FontSize',8)
    ylim([lb,ub])
    box off

    subplot(1,2,1)
    hold on
    legend('Animal recording','Model response','FontSize',8)
    ylabel('AP/s','FontSize',8)
    xlabel('time (s)','FontSize',8)
    title('(a) Sensory discharge','FontSize',8)
    box off
    
    hOneGroupFig.Position(3) = 700;
    hOneGroupFig.Position(4) = 225;
    
    set(hOneGroupFig,'renderer','Painters')
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TWO-GROUP INVERSION METHOD
% FIGURE 10
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

useBoth = true;
 
for i=5:8

    hTwoGroupFig = figure;

    lb = -3;
    ub = 3;
    p0 = .1;
    xTol = 1e-4;
    fTol = 1e-10;
    kMax = 20;

    t = cleanedData{i}.t(1:end-1);
    u6A = cleanedData{i}.u(2:end);
    u6B = -cleanedData{i}.u(2:end);
    
    sens6A = cleanedData{i}.sens6A(2:end);
    sens6B = cleanedData{i}.sens6B(2:end);
    
    nSamps = length(t);
    uCalc = zeros(nSamps,2); %column 1 is lower limit, column 2 is upper limit
    
    if useBoth
        
        tauComb = [tau6A,tau6B];
        aComb = [a6A,a6B];
        bComb = [b6A,b6B];
        cComb = [c6A,c6B];
        dComb = [d6A,d6B];
        sensComb = [sens6A,sens6B];
        
        x06A = u6A(1) - (sens6A(1) - c6A*u6A(1) - d6A)/a6A;

        x06B = u6B(1) - (sens6B(1) - c6B*u6B(1) - d6B)/a6B;
        
        x0CombUpper = [x06A,x06B];
        x0CombLower = [x06A,x06B];
        
        toInvert = [false,true];
    else
        %6B only
        x0CombUpper = x02;
        x0CombLower = x02;
        
        tauComb = xf6B(1)*1e-3;
        aComb = xf6B(2)*1e3;
        bComb = xf6B(3);
        cComb = xf6B(4);
        dComb = xf6B(5);
        sensComb = sens6B;
        
        toInvert = true;
    end
    toPlot = false;

    for j=2:nSamps
        fprintf('Both, 6A = %f and 6B = %f\n',sens6A(j),sens6B(j))

        uLower = [uCalc(:,1),uCalc(:,1)];
        uUpper = [uCalc(:,2),uCalc(:,2)];

        fvalTol = .5; 

        if j == 2
            if useBoth
                %first loop, use initial value
                x0CombUpper = [x06A,x06B];
                x0CombLower = [x06A,x06B];
            else
                %first loop, use initial value
                x0CombUpper = x06B;
                x0CombLower = x06B;
            end
        else
            %not first loop, use previous final value
            if useBoth
                x0CombUpper = [xHi6A(end),xHi6B(end)];
                x0CombLower = [xLo6A(end),xLo6B(end)];
            else
                %6B only
                x0CombUpper = xHi6B(end);
                x0CombLower = xLo6B(end);
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        gCombLower = @(x)loadMatchingFuncMulti(x,j,tauComb,aComb,bComb,cComb,dComb,t,uLower,sensComb,x0CombLower,toInvert,toPlot);
        x0 = uCalc(j-1,1);

        if gCombLower(x0) == 0
            %do nothing
        else
            try x0 = fzero(gCombLower,[lb,ub]);
            catch
                k = 0;
                x0Vec = NaN(kMax+1,1);
                fVec = NaN(kMax+1,1);
                badFit = true;
                while badFit && (k < 20)
                    [x0,fvalLower,exitFlagLower] = fminsearch(gCombLower,x0+(-1)^k*k*.025);
                    
                    x0Vec(k+1) = x0;
                    fVec(k+1) = fvalLower;
                    
                    fprintf('Lower initial x: %f\n',x0+(-1)^k*k*.025)
                    fprintf('Lower exit value: %f\n',fvalLower)

                    k = k + 1;
                    if abs(fvalLower) < fvalTol
                        badFit = false;
                    else
                        fprintf('STILL BAD FIT.\n')
                        gCombLower = @(x)loadMatchingFuncMulti(x,j,tauComb,aComb,bComb,cComb,dComb,t,uLower,sensComb,x0CombLower,toInvert,toPlot);
                    end
                end
                if badFit
                    %Find the x0 value that had the lowest assoc. f value
                    [~,minInd] = min(fVec);
                    x0 = x0Vec(minInd);
                end
            end
        end

        uCalc(j,1) = x0;

        %%%%%%%%%%%%%%%%%%%%%%%%%%%


        gCombUpper = @(x)loadMatchingFuncMulti(x,j,tauComb,aComb,bComb,cComb,dComb,t,uUpper,sensComb,x0CombUpper,toInvert,toPlot);
        x0 = uCalc(j-1,2);

        if gCombUpper(x0) == 0
            %do nothing
        else
            try x0 = fzero(gCombUpper,[lb,ub]);
            catch
                k = 0;
                x0Vec = NaN(kMax+1,1);
                fVec = NaN(kMax+1,1);
                badFit = true;
                while badFit && (k < 20)
                    [x0,fvalUpper,exitFlagUpper] = fminsearch(gCombUpper,x0+(-1)^k*k*.025);
                    
                    x0Vec(k+1) = x0;
                    fVec(k+1) = fvalUpper;
                    
                    fprintf('Upper initial x: %f\n',x0+(-1)^k*k*.025)
                    fprintf('Upper exit value: %f\n',fvalUpper)

                    k = k + 1;
                    if abs(fvalUpper) < fvalTol
                        badFit = false;
                    else
                        fprintf('STILL BAD FIT.\n')
                        gCombUpper = @(x)loadMatchingFuncMulti(x,j,tauComb,aComb,bComb,cComb,dComb,t,uUpper,sensComb,x0CombUpper,toInvert,toPlot);
                    end
                end
                if badFit
                    %Find the x0 value that had the lowest assoc. f value
                    [~,minInd] = min(fvalUpper);
                    x0 = x0Vec(minInd);
                end
            end
        end

        uCalc(j,2) = x0;
    
        [yLo6A,xLo6A] = simulateMinusLowpassPL(uCalc(1:j,1),tau6A,t(1:j),a6A,b6A,c6A,d6A,x06A);
        [yHi6A,xHi6A] = simulateMinusLowpassPL(uCalc(1:j,2),tau6A,t(1:j),a6A,b6A,c6A,d6A,x06A);
        [yLo6B,xLo6B] = simulateMinusLowpassPL(-uCalc(1:j,1),tau6B,t(1:j),a6B,b6B,c6B,d6B,x06B);
        [yHi6B,xHi6B] = simulateMinusLowpassPL(-uCalc(1:j,2),tau6B,t(1:j),a6B,b6B,c6B,d6B,x06B);
        x0CombUpper = [xHi6A(end),xHi6B(end)];
        x0CombLower = [xLo6A(end),xLo6B(end)];

        figure(hTwoGroupFig)
        subplot(2,2,1)
        cla
        plot(t(1:j),sens6A(1:j),'linewidth',2)
        hold on
        plot(t(1:j),yLo6A,':','linewidth',2)

        %6B responses
        subplot(2,2,3)
        cla
        plot(t(1:j),sens6B(1:j),'linewidth',2)
        hold on
        plot(t(1:j),yLo6B,':','linewidth',2)
        
        forceAx = subplot(2,2,[2,4]);
        cla
        plot(t(1:j),-u6A(1:j),'linewidth',2)
        hold on
        plot(t(1:j),-uCalc((1:j),1),'--','linewidth',2)
        plot(t(1:j),-uCalc((1:j),2),':','linewidth',2)
        drawnow
    end

    lineColors = lines(3);

    figure(hTwoGroupFig)
    subplot(forceAx);
    hold on
    ylim([-.5,3])
    legend('Actual force','Lower bound force','Upper bound force','location','northwest','Fontsize',8)
    title('(c) Estimated and actual stimulus force','FontSize',8)
    xlabel('Time (s)','Fontsize',8)
    ylabel('mN','FontSize',8)
    box off

    subplot(2,2,1)
    hold on
    ylim([0,250])
    legend('Animal recording','Model response','location','south','FontSize',8)
    title('(a) 6A discharge','Fontsize',8)
    ylabel('AP/s','FontSize',8)
    box off

    subplot(2,2,3);
    hold on
    ylim([0,250])
    legend('Animal recording','Model response','location','south','FontSize',8)    
    title('(b) 6B discharge','Fontsize',8)
    ylabel('AP/s','FontSize',8)
    xlabel('Time (s)','FontSize',8)
    box off
    
    hTwoGroupFig.Position(3) = 700;
    hTwoGroupFig.Position(4) = 225;
    
    forceAx.Position(2) = .1725;
    forceAx.Position(4) = .7526;
    
    set(hTwoGroupFig,'renderer','Painters')

end

