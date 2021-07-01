% tuneCockroach
% Nicholas Szczecinski
% Department of Mechanical and Aerospace Engineering
% West Virginia University
% 2 April 2021

clear
close all
clc

%Stimulus waveform used in Zill et al. 2018, Figure 9
ctr = importdata('CTrTorqueWaveform.csv');

%Ramp-And-Hold (RAH) stimulus
RAHanimal = load('rampAndHold.mat');
RAHanimal.heightsMapped(isnan(RAHanimal.heightsMapped)) = [];
RAHanimal.heightsMapped = [0;RAHanimal.heightsMapped];
n = sum(~isnan(RAHanimal.heightsMapped));
tRAH = linspace(0,100,n);

TWanimal = load('torqueWaveform.mat');
TWanimal.heightsMapped(isnan(TWanimal.heightsMapped)) = [];
TWanimal.heightsMapped = [0;TWanimal.heightsMapped];
tTW = linspace(0,100,sum(~isnan(TWanimal.heightsMapped)));

t = ctr.data(:,1);
dt = mean(diff(t));
numSteps = length(t);
tmax = max(t);

%Initialize parameter values
tau = 4e-3;
a = 700;
b = 2.3;
c = 70; 
d = -60; 

%U1: Ramp and hold
u1 = zeros(numSteps,1);
A = 1.2;
T = .17*tmax;
tStart = 0;
tEnd = .62;
u1(t >= tStart) = min(A,A/T*t(1:sum(t >= tStart)));
u1(t >= tEnd) = max(0,A-A/T*t(1:sum(t >= tEnd)));

h1 = 1e-4;
t1 = (0:h1:max(t))';
u11 = interp1(t,u1,t1);
[y11,x11,du11dt] = simulateMinusLowpassPL(u11,tau,t1,a,b,c,d,0);

x0 = [50;.5];

lb = [1;0;1;0;-100];
ub = [100;2;5;100;100];

t1t = linspace(0,tmax,length(RAHanimal.heightsMapped))';
u1t = interp1(t,u1,t1t);
h1t = mean(diff(t1t));

g = @(x) max(0,x(1)*fgl_deriv( x(2), u1t, h1t ));
f = @(x) sum( (g(x) - RAHanimal.heightsMapped).^2);

figure
hold on
plot(t1t,RAHanimal.heightsMapped)
plot(t1t,50*fgl_deriv( .5, u1t, h1t ),'--')

options = optimoptions('fmincon','OptimalityTolerance',1e-12,'StepTolerance',1e-12,'UseParallel',false,'Display','iter');
xf = fmincon( f,x0,[],[],[],[],lb,ub,[],options);

% Values used in the manuscript:
% xf = [
% 
%     3.8595;
%     0.7076;
%     2.2615;
%    54.2870;
%   -41.1551];

a = xf(1);
b = xf(2);

cockroachParamsFrac.a = a;
cockroachParamsFrac.b = b;

save('cockroachParamsFrac.mat','cockroachParamsFrac');

figure
hold on
plot(t1t,RAHanimal.heightsMapped)
plot(t1t,max(0,a*fgl_deriv( b, u1t, h1t )),'--')


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOT FIGURE 2 FROM THE MANUSCRIPT.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%U2: Torque waveform
%Factor of 0.6 scales the stimulus as presented in Zill et al. 2018, Fig 9.
u2 = 0.6*ctr.data(:,2);

% [y1,x1,du1dt] = simulateMinusLowpassPL(u1,tau,t,a,b,c,d,0);
% [y2,x2,du2dt] = simulateMinusLowpassPL(u2,tau,t,a,b,c,d,lpf0(xf));

y1 = max(0,a*fgl_deriv( b, u1, mean(diff(t)) ));
y2 = max(0,a*fgl_deriv( b, u2, mean(diff(t)) ));

h = figure;
h.Position(2) = 100;
h.Position(3) = 750;
h.Position(4) = 450;

colors = lines(6);

meanAbsErrFrac = NaN(2,1);

spa = subplot(2,4,1);
p3 = plot(t/tmax*100,100*u1,'k','linewidth',2);
hold on
p1 = plot(tRAH,RAHanimal.heightsMapped,'linewidth',2,'color',colors(1,:));
p2 = plot(t/tmax*100,y1,':','linewidth',2,'color',colors(2,:));
title({'(a) Response to ramp-','and-hold stimulus'},'Fontsize',8)
ylim([0,200])
xlabel('Percent of stimulus','FontSize',8)
xticks(0:20:100)
box off
lgd = legend( [p1,p2,p3],{'Animal response (Hz)','Model response (Hz)','Force (100x mN)'},'Location','Southoutside');
legend('boxoff')

meanAbsErrFrac(1) = mean(abs(RAHanimal.heightsMapped - interp1(t/tmax*100,y1,tRAH')));

spb = subplot(2,4,2);
p6 = plot(t/tmax*100,100*u2,'k','linewidth',2);
hold on
p5 = plot(tTW,TWanimal.heightsMapped,'linewidth',2,'color',colors(1,:));
p4 = plot(t/tmax*100,y2,':','linewidth',2,'color',colors(2,:));

title({'(b) Response to','naturalistic stimulus'},'Fontsize',8)
ylim([0,200])
yticks([])
xlabel('Percent of stimulus','FontSize',8)
xticks(0:20:100)
box off
legend( [p4,p5,p6],{'Animal Response (Hz)','Model Response (Hz)','Force (100x mN)'},'Location','Southoutside')
legend('hide')

meanAbsErrFrac(2) = mean(abs(TWanimal.heightsMapped - interp1(t/tmax*100,y2,tTW')));

barColors = lines(4);
barColors(1:2,:) = [];
spc = subplot(2,4,[3,4]);
bar(t/tmax*100,y1,.5,'LineStyle','none','FaceColor',barColors(1,:))
hold on
ttc = title('(c) Comparison between sensory responses','Color','k','FontSize',8);
ttc.Position(2) = 212;

bar(t/tmax*100+.5,y2,.5,'LineStyle','none','FaceColor',barColors(2,:))
ylim([0,200])
xlim([0,100])
xlabel('Percent of stimulus','FontSize',8)
ylabel('Afferent response (Hz)','FontSize',8)
xticks(0:20:100)

legend(['Model ramp-',newline,'and-hold response'],['Model naturalistic',newline,'response'],'Location','Northeast')
legend('boxoff')

save('cockroachFitErrFrac.mat','meanAbsErrFrac');

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
% MAKE SURE THIS TUNING GENERALLY WORKS BY APPLYING MORE RAMP STIMULI
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Load data from cockroach proximal tibial CS
diffRates = load('sameAmplitudeDifferentRate.csv');
Tvec = diffRates(1,2:end);
A = 1.66;

Tramp = NaN(4,1);
Aramp = NaN(4,1);
tRamp = NaN(21,4);

for i=1:4

    Tramp(i) = Tvec(i);
    Aramp(i) = A;
    dtRamp = Tvec(i)/20;
    tRamp(:,i) = (0:dtRamp:Tvec(i))';
    
    if i == 3
        Tramp(3) = .225;
    elseif i == 4
        Tramp(4) = .125;
    end
    
    uRamp = A/Tvec(i)*tRamp;
    response(:,i) = diffRates(:,1+i); %#ok<SAGROW>
    response(1,i) = 0; %#ok<SAGROW>
end

%Simulate model response to each ramp and plot it to produce figure 2.
subplotIndices = [5,6,7,8];

meanAbsErrRampsFrac = NaN(4,1);

for i=1:4    
    
    dt = 1e-2;
    tSim = (0:dt:1)';
    
    u = Aramp(i)/Tramp(i)*tSim;
    u(u > Aramp(i)) = Aramp(i);
    
    %get the right time and u
%     y = simulateMinusLowpassPL(u,tau,tSim,a,b,c,d,0);
    y = max(0,a*fgl_deriv( b, u, mean(diff(tSim)) ));
    
    figure(h)
    hold on
    subplot(2,4,subplotIndices(i))
    hold on
    p3 = plot(tSim,100*u,'k','linewidth',2);
    p1 = stairs(tRamp(:,i),response(:,i),'linewidth',2,'color',colors(1,:));
    p2 = plot(tSim,y,':','linewidth',2,'color',colors(2,:));
    
    ylim([0,300])
    
    if i == 1
        tt = title('(d) Response to additional ramp-and-hold stimuli','fontsize',8,'HorizontalAlignment','left');
        tt.Position(1) = -.2;
        tt.Position(2) = 312;
    end
    xlabel('t (s)','FontSize',8)
    box off
    
    
    if i == 1 
        ylabel('y (Hz)')
    else
        yticks([])
    end
end

spc.Position(3) = .3;
spc.Position(1) = .6;

hY = spa.Position(4);
lgd.Position = [0.1250 0.5210 0.1993 0.0965];
spa.Position = [0.1300 0.6940 0.1521 0.2239];
spb.Position = [0.3361 0.6940 0.1521 0.2239];

save('cockroachFitErrorRampsFrac.mat','meanAbsErrRampsFrac')