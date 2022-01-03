% characterizeCockroach
% Nicholas Szczecinski
% Department of Mechanical and Aerospace Engineering
% West Virginia University
% 2 April 2021

clear
close all
clc

dt = 1e-4;
t = (0:dt:3)';

numSteps = length(t);
tmax = max(t);

load('cockroachParams');
tau = cockroachParams.tau;
a = cockroachParams.a;
b = cockroachParams.b;
c = cockroachParams.c;
d = cockroachParams.d;

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% QUANTIFY HYSTERESIS.
% PLOT FIGURE 5 FROM THE MANUSCRIPT.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

th = (0:dt:27)';
uh = zeros(size(th));
duh = zeros(size(th));
Th = 1.5;
T3 = 0.2246927273;
Ah = 0.5;
Acurrent = 0;
Amax = 4;
rampUp = false;
rampDn = false;
ascending = true;
iStartTonic = NaN(size(th));
iEndTonic = NaN(size(th));
j = 1;
for i=2:length(th)
    if mod(th(i),Th) == 0
        
        iEndTonic(j) = i;
        
        j = j + 1;
        if ascending
            Acurrent = Acurrent + Ah;
            %start a new ramp
            rampUp = true;
            rampDn = false;
            if Acurrent > Amax
                Acurrent = Amax;
                rampUp = false;
                rampDn = false;
                ascending = false;
            end
        else
            Acurrent = Acurrent - Ah;
            rampDn = true;
            rampUp = false;
        end
    end
    
    if rampUp
        uh(i) = Acurrent - Ah + Ah/T3*mod(th(i),Th);
        duh(i) = Ah/T3;
        if uh(i) > Acurrent
            iStartTonic(j) = i;
            rampUp = false;
        end
    elseif rampDn
        uh(i) = Acurrent + Ah - Ah/T3*mod(th(i),Th);
        duh(i) = -Ah/T3;
        if uh(i) < Acurrent 
            iStartTonic(j) = i;
            rampDn = false;
        end
    else
        uh(i) = Acurrent;
    end
end

dsFactor = 100; %downsample by a factor of dsFactor

%Plot staircase stimulus
g = figure;
sp3 = subplot(6,2,[1,3]);
plot(downsample(th,dsFactor),downsample(uh,dsFactor),'k','linewidth',1)
ylabel('u (mN)','fontsize',8)
xlim([0,max(th)])
title({'(a) ''Staircase'' stimulus'},'Fontsize',8)
xticklabels([])
box off
drawnow

%Simulate response to the staircase
[yh,xh,dudt] = simulateMinusLowpassPL(uh,tau,th,a,b,c,d,0);
hold on

%Plot the response to the staircase in plain black (individual 'steps' will
%be added below).
sp1 = subplot(6,2,[5,7]);
plot(downsample(th,dsFactor),downsample(yh,dsFactor),'k','linewidth',.5)
hold on
xlim([0,max(th)])
ylabel('y (Hz)','fontsize',8);
title({'(b) Model response to ''staircase'' stimulus'},'Fontsize',8)
xticklabels([])
box off

%Calculate the response of the 'dynamics free' model, wherein y explicitly
%relies on du/dt.
yNoDyn = max(0,a/c*sign(dudt).*abs(dudt).^b + c*uh + d);
spNoDyn = subplot(6,2,[9,11]);
plot(downsample(th,dsFactor),downsample(yNoDyn,dsFactor),'k','linewidth',.5)
hold on
xlim([0,max(th)])
title('(c) Response of dynamics-free model','Fontsize',8)
xlabel('time (s)','fontsize',8)
ylabel('y (Hz)','fontsize',8);
box off

%Define a color scheme for plotting the response to each step.
colorGrad = parula(j-1);

avgAct = NaN(1,j-1);
avgForce = NaN(1,j-1);
avgNoDynAct = NaN(1,j-1);

for i=1:(j-1)
    if iEndTonic(i) > iStartTonic(i) && iStartTonic(i) > 0
        inds = (iStartTonic(i)+1):(iEndTonic(i)-1);
        avgAct(i) = mean(yh(inds));
        avgForce(i) = mean(uh(inds));
        
        %Plot the time following each step in color, over the plain black
        %response plot
        subplot(sp1)
        plot(downsample(th(inds),dsFactor),downsample(yh(inds),dsFactor),'color',colorGrad(i,:),'linewidth',2);
        
        avgNoDynAct(i) = mean(yNoDyn(inds));
        
        %Plot the time following each step in color, over the plain black
        %response plot
        subplot(spNoDyn)
        plot(downsample(th(inds),dsFactor),downsample(yNoDyn(inds),dsFactor),'color',colorGrad(i,:),'linewidth',2);
    end
end

%To simplify summary plots in Fig 5 D and E, remove points during which the
%average force was 0.
colorGrad(avgForce == 0,:) = [];
avgAct(avgForce == 0) = [];
avgNoDynAct(avgForce == 0) = [];
avgForce(avgForce == 0) = [];

%Boolean mask of steps for which the force increased to that level (as
%opposed to decreased).
increasing = [true,(diff(avgForce) > 0)];

%Plot summary of responses to reveal hysteresis
sph = subplot(6,2,2:2:6);
scatter(avgForce(increasing),avgAct(increasing),[],colorGrad(increasing,:),'^','filled')
hold on
scatter(avgForce(~increasing),avgAct(~increasing),[],colorGrad(~increasing,:),'v','filled')
title({'(d) Model response exhibits hysteresis upon repeated loading'},'Fontsize',8,'Position',[2,203,0])
ax = gca;
ax.YLim(1) = 0;
ax.XLim(1) = 0;
xticklabels([])
ylabel('y (Hz)','fontsize',8)
box on
ax.Position = [.57,.5568,.3347,.3344];

%Plot summary of responses to reveal hysteresis
spNoDynSummary = subplot(6,2,8:2:12);
scatter(avgForce(increasing),avgNoDynAct(increasing),[],colorGrad(increasing,:),'^','filled')
hold on
scatter(avgForce(~increasing),avgNoDynAct(~increasing),[],colorGrad(~increasing,:),'v','filled')
title({'(e) Without adaptation, hysteresis is not present'},'Fontsize',8,'Position',[2,202,0])
ax = gca;
ax.YLim(1) = 0;
ax.XLim(1) = 0;
ax.Position = [.57,.1063,.3347,.348]; 
box on

xlabel('u (mN)','fontsize',8)
ylabel('y (Hz)','fontsize',8)

sgtitle('Model Response to "Staircase" Stimulus Exhibits Hysteresis','Fontsize',10)
g.Position = [394.6000 157.8000 700.0000 500.0000];

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CHARACTERIZE SENSITIVITY TO STIMULUS AMPLITUDE AND OFFSET
% PLOT FIGURE 3 FROM THE MANUSCRIPT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

hStepAndRate = figure;

amplitudes = [0.6;1.3;1.9;2.8];
numAmplitudes = length(amplitudes);
durations = [0.7;0.6;0.5;0.4;0.3;0.2;0.1;0.05;0.02];
numDurations = length(durations);
rate = NaN(numAmplitudes,numDurations);
maxFreq = NaN(numAmplitudes,numDurations);
meanFreq = NaN(numAmplitudes,numDurations);
lineColors = lines(4);
ratesExamplePlotIndex = 6;
for i=1:numAmplitudes
    A = amplitudes(i);
    for j=1:numDurations
        T = durations(j);
        tOff = t;
        uRelease = min(A/T*tOff,A);
        
        rate(i,j) = A/T;
        
        [y,x] = simulateMinusLowpassPL(uRelease,tau,tOff,a,b,c,d,0);
        maxFreq(i,j) = max(y);
        meanFreq(i,j) = mean(y(tOff >= T));
        
        if j == ratesExamplePlotIndex
            subplot(4,2,2*i-1)
            plot(tOff,100*uRelease,'k','linewidth',1)
            hold on
            plot(tOff,y,'color',lineColors(i,:),'linewidth',1)
            ylim([0,400])
            xlim([0,1.5])
            box off
            if i < 4
                xticks([])
            else
                xlabel('time (sec)','FontSize',8)
            end
        end
    end
end

mkrCell = {'s','o','^','d'};
fillVec = [false,true,false,true];


xmin = 0.5;
xmax = 40;
xlims = [xmin;xmax];
sp = NaN(4,1);
slopes = NaN(4,1);
legLabels = cell(1,4);

subplot(4,2,2:2:8)
hold on
for i=1:numAmplitudes
    %Only keep trials where the rate is between 0.5 and 40 mN/sec. Otherwise,
    %ignore them.
    bm = (rate(i,:) >= xmin & rate(i,:) <= xmax);
    
    B = [log10(rate(i,bm)'),ones(size(rate(i,bm)'))];
    C = log10(maxFreq(i,bm)');
    fitCoeffs = (B'*B)^-1*B'*C;
    slopes(i) = fitCoeffs(1);
    
    testDomain = [log10(xlims),[1;1]];
    
    plot(10.^testDomain(:,1),10.^(testDomain*fitCoeffs),'color',lineColors(i,:));
    hold on
    drawnow
    
    if fillVec(i)
        scatter(rate(i,ratesExamplePlotIndex),maxFreq(i,ratesExamplePlotIndex),72,mkrCell{i},'k')
        hold on
        sp(i) = scatter(rate(i,bm),maxFreq(i,bm),mkrCell{i},'filled','MarkerFaceColor',lineColors(i,:));
    else
        scatter(rate(i,ratesExamplePlotIndex),maxFreq(i,ratesExamplePlotIndex),72,mkrCell{i},'k')
        hold on
        sp(i) = scatter(rate(i,bm),maxFreq(i,bm),mkrCell{i},'filled','MarkerFaceColor','w','MarkerEdgeColor',lineColors(i,:));
    end
    legLabels{i} = sprintf('%1.1f mN amplitude (k = %0.2f)',amplitudes(i),slopes(i));
end
ax = gca;
ax.XScale = 'log';
ax.YScale = 'log';
ylim([10,1000])
xlim(xlims)
xlabel('Rate of change of force, du/dt (mN/sec)','FontSize',8)
ylabel('Maximum discharge (Hz)','FontSize',8)
hStepAndRate.Position(3) = 700;
hStepAndRate.Position(4) = 300;
ax.Position(2) = .16;
ax.Position(4) = .69;
box on

legend(sp,legLabels,'location','southeast','fontsize',8)

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RESPONSE TO RAMPS WITH VARIOUS OFFSETS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

hOffset2 = figure;

offsets = -[0;-0.2;-0.4;-0.6];
offsets = [0.1;0.7;1.3;1.5];
A = 1;
numOffsets = length(offsets);
durations = [0.7;0.6;0.5;0.4;0.3;0.2;0.1;0.05;0.02];
numDurations = length(durations);
rate = NaN(numOffsets,numDurations);
maxFreq = NaN(numOffsets,numDurations);
meanFreq = NaN(numOffsets,numDurations);
lineColors = lines(4);
tOnset = 0.5;
offsetRateExamplePlotIndex = 6;
for i=1:numOffsets
    B = abs(offsets(i));
    sB = sign(offsets(i));
    for j=1:numDurations
        T = durations(j);
       
        uRelease = sB*B + min(max(0,A/T*(t-tOnset)),A);
        
        rate(i,j) = A/T;
        
        [y,x] = simulateMinusLowpassPL(uRelease,tau,t,a,b,c,d,sB*B);
        maxFreq(i,j) = max(y);
        meanFreq(i,j) = mean(y(t >= tOnset+T));
        
        if j == offsetRateExamplePlotIndex
            subplot(4,2,2*i-1)
            plot(t,100*uRelease,'k','linewidth',1)
            hold on
            plot(t,y,'color',lineColors(i,:),'linewidth',1)
            ylim([0,250])
            xlim([0,1.5])
            box off
            if i < 4
                xticks([])
            else
                xlabel('time (sec)','FontSize',8)
            end
        end
    end
end

mkrCell = {'s','o','^','d'};
fillVec = [false,true,false,true];


xmin = 0.5;
xmax = 40;
xlims = [xmin;xmax];
sp = NaN(4,1);
slopes = NaN(4,1);
legLabels = cell(1,4);

subplot(4,2,2:2:8)
hold on
for i=1:numOffsets
    %Only keep trials where the rate is between 0.5 and 40 mN/sec. Otherwise,
    %ignore them.
    bm = (rate(i,:) >= xmin & rate(i,:) <= xmax);
    
    B = [log10(rate(i,bm)'),ones(size(rate(i,bm)'))];
    C = log10(maxFreq(i,bm)');
    fitCoeffs = (B'*B)^-1*B'*C;
    slopes(i) = fitCoeffs(1);
    
    testDomain = [log10(xlims),[1;1]];
    
    plot(10.^testDomain(:,1),10.^(testDomain*fitCoeffs),'color',lineColors(i,:));
    hold on
    drawnow
    
    if fillVec(i)
        scatter(rate(i,offsetRateExamplePlotIndex),maxFreq(i,offsetRateExamplePlotIndex),72,mkrCell{i},'k')
        hold on
        sp(i) = scatter(rate(i,bm),maxFreq(i,bm),mkrCell{i},'filled','MarkerFaceColor',lineColors(i,:));
    else
        scatter(rate(i,offsetRateExamplePlotIndex),maxFreq(i,offsetRateExamplePlotIndex),72,mkrCell{i},'k')
        hold on
        sp(i) = scatter(rate(i,bm),maxFreq(i,bm),mkrCell{i},'filled','MarkerFaceColor','w','MarkerEdgeColor',lineColors(i,:));
    end
    legLabels{i} = sprintf('%1.1f mN offset (k = %0.2f)',offsets(i),slopes(i));
end
ax = gca;
ax.XScale = 'log';
ax.YScale = 'log';
ylim([10,1000])
xlim(xlims)
xlabel('Rate of change of force, du/dt (mN/sec)','FontSize',8)
ylabel('Maximum discharge (Hz)','FontSize',8)
hOffset2.Position(3) = 700;
hOffset2.Position(4) = 300;
ax.Position(2) = .16;
ax.Position(4) = .69;
box on

legend(sp,legLabels,'location','southeast','fontsize',8)

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DISCHARGE DECAYS OVER TIME ACCORDING TO A POWERLAW
% PLOT FIGURE 4
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

hDecay = figure;
tDecay = [0;logspace(-4,2,1000)'];
A = .1*[1;2;3;4];

for i=1:4
    uRelease = A(i)+zeros(size(tDecay));

    [y,x] = simulateMinusLowpassPL(uRelease,tau,tDecay,a,b,0,0,0);

    beta = 1-b;
    s = 1/beta;
    B = (-beta)^(1/beta);
    delta_t = tau*(A(i)/B)^(1-b);
    xh = B*((tDecay+delta_t)/tau).^s;

    yExact = a*(xh); % + c*A(i) + d;

    subplot(4,2,2*i-1)
    plot(tDecay+delta_t,y,'k','linewidth',1);
    hold on
    plot(tDecay+delta_t,yExact,':','linewidth',1+(i-1)/2,'color',lineColors(i,:))
    xlim([0,2])
    box off
    if i == 4
        xlabel('time (s)','fontsize',8)
    else
        xticks([])
    end
    
    ax = subplot(4,2,2:2:8);
    plTemp = plot(tDecay+delta_t,y,'k','linewidth',1);
    hold on
    if i == 1
        plVec(1) = plTemp;
        plVec(2) = plot([1,2],[1,2],'linestyle','none');
    end
    plVec(i+2) = plot(tDecay+delta_t,yExact,':','linewidth',1+(i-1)/2,'color',lineColors(i,:)); %#ok<SAGROW>
    xlim([1e-2,max(tDecay+delta_t)])
    if i == 4
        xlabel('time (s)','fontsize',8)
        ylabel('discharge rate, y (Hz)','fontsize',8)
    end
    box on
end
ax.XScale = 'log';
ax.YScale = 'log';
ax.Position(2) = .16;
ax.Position(4) = .69;
legend(plVec,{'simulated responses','analytical responses:','A = 0.1 mN','A = 0.2 mN','A = 0.3 mN','A = 0.4 mN'},'fontsize',8,'location','southwest')

hDecay.Position(3) = 700;
hDecay.Position(4) = 300;




