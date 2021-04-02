% invertCockroach
% Nicholas Szczecinski
% Department of Mechanical and Aerospace Engineering
% West Virginia University
% 2 April 2021

clear
close all
clc

load('cockroachParams');
tau = cockroachParams.tau;
a = cockroachParams.a;
b = cockroachParams.b;
c = cockroachParams.c;
d = cockroachParams.d;

ctr = importdata('CTrTorqueWaveform.csv');

t = ctr.data(:,1);
dt = mean(diff(t));
numSteps = length(t);
tmax = max(t);

%When testing the inversion algorithm, it was sometimes useful to ask the
%model to invert output from the forward model. The forward model's output
%is less noisy than the animal recordings, ensuring that a perfect solution
%can be found. However, we also show in the manuscript that this approach
%works with noisy animal data. Thus, 'useExactAnimalRecordings' is set to
%true, but you may change the 'desired' output by changing it to false.
useExactAnimalRecordings = true;

%U1: Ramp and hold
u1 = zeros(numSteps,1);
A = 1.2;
T = .17*tmax;
tStart = 0;
tEnd = .62;
u1(t >= tStart) = min(A,A/T*t(1:sum(t >= tStart)));
u1(t >= tEnd) = max(0,A-A/T*t(1:sum(t >= tEnd)));

%U2: Torque waveform
u2 = 0.6*ctr.data(:,2);

[y1,x1,du1dt] = simulateMinusLowpassPL(u1,tau,t,a,b,c,d,0);
[y2,x2,du2dt] = simulateMinusLowpassPL(u2,tau,t,a,b,c,d,0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RAMP AND HOLD, FIGURE 8 A and B
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if useExactAnimalRecordings
    m = 1;
    RAHanimal = load('rampAndHold.mat');
    n = sum(~isnan(RAHanimal.heightsMapped));
    tRAH = linspace(0,tmax,n+1)';
    RAHanimal.heightsMapped(isnan(RAHanimal.heightsMapped)) = [];
    yRAH = [0;RAHanimal.heightsMapped];
    tSamp = linspace(0,tmax,floor((n+1)/m))';
else
    m = 1; %#ok<UNRCH>
    n = sum(~isnan(y1));
    tRAH = linspace(0,tmax,n+1)';
    y1(isnan(y1)) = [];
    yRAH = [0;y1];
    tSamp = linspace(0,tmax,floor((n+1)/m))';
    uDes = [0;u1];
end

yDesRAH = yRAH;

nSamps = length(tSamp);
uSamp = zeros(nSamps,2);

x01 = 0;
x02 = 0;

%Parameters for the cornerSearch. p0 is the initial step length of the
%algorithm. xTol establishes the minimum step length allowed (if reached,
%the algorithm halts). fTol establishes how much of a deviation from 0 is
%considered a 'corner'.
p0 = .1;
xTol = 1e-4;
fTol = 1e-10;

%Lower bounds and upper bounds for the cornerSearch.
lb = -eps;
ub = 3;

h1 = figure;

for i=2:nSamps
    
    x0 = uSamp(i-1,1);
    imax = i;
    tLoop = tRAH(1:imax);
    yDesRAH = yRAH(1:imax);

    g1 = @(x)loadMatchingFunc(x,i,tau,a,b,c,d,tSamp,uSamp(:,1),yDesRAH,x01);
    
    if g1(x0) == 0
        %do nothing
    else
        x0 = fzero(g1,[lb,ub]);
    end
    
    [uSamp(i,1),~,~,converged] = cornerSearch(g1,x0,p0,lb,ub,100,xTol,fTol,false);
    if ~converged
        error('cornerSearch did not converge. It likely reached the maximum number of iterations. Consider a larger value for p0.')
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    x0 = uSamp(i-1,2);

    g2 = @(x)loadMatchingFunc(x,i,tau,a,b,c,d,tSamp,uSamp(:,2),yDesRAH,x02);
    
    if g2(x0) == 0
        %do nothing
    else
        x0 = fzero(g2,[lb,ub]);
    end
    
    [uSamp(i,2),~,~,converged] = cornerSearch(g2,x0,p0,lb,ub,100,xTol,fTol,true);
    if ~converged
        error('cornerSearch did not converge. It likely reached the maximum number of iterations. Consider a larger value for p0.')
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    figure(h1)
    clf
    subplot(1,2,1)
    plot(tRAH(1:i),yDesRAH,'linewidth',2)
    hold on
    [yLo,xLo] = simulateMinusLowpassPL(uSamp(1:i,1),tau,tSamp(1:i),a,b,c,d,0);
    plot(tLoop,yLo,':','linewidth',2)
    [yHi,xHi] = simulateMinusLowpassPL(uSamp(1:i,2),tau,tSamp(1:i),a,b,c,d,0);
    
    x01 = xLo(end);
    x02 = xHi(end);
    
    subplot(1,2,2)
    plot(t,u1,'linewidth',2)
    hold on
    plot(tSamp(1:i),uSamp(1:i,1),'--','linewidth',2)
    plot(tSamp(1:i),uSamp(1:i,2),':','linewidth',2)
    drawnow
    
end

uInterped = @(u) pchip(tSamp,u,tRAH);

figure(h1)
subplot(1,2,2)
hold on
legend('actual force','lower bound force','upper bound force','location','south','FontSize',8)
title('Applied force','FontSize',8)
xlabel('time (s)','FontSize',8)
ylabel('mN','FontSize',8)
ylim([0,1.5])
box off

subplot(1,2,1)
hold on
title('Sensory discharge','FontSize',8)
ylabel('Hz','FontSize',8)
xlabel('time (s)','FontSize',8)
ylim([0,200])
box off

if useExactAnimalRecordings
    sgtitle('(a) Estimated applied force from discharge pattern (animal)','Fontsize',10)
else
    sgtitle('(a) Estimated applied force from discharge pattern (simulation)','Fontsize',10) %#ok<UNRCH>
end
h1.Position(3) = 700;
h1.Position(4) = 200;
set(h1,'renderer','Painters')

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TORQUE WAVEFORM/NATURALISTIC STIMULUS, FIGURE 8C and D 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if useExactAnimalRecordings
    m = 1;
    TWanimal = load('torqueWaveform.mat');
    n = sum(~isnan(TWanimal.heightsMapped));
    tTW = linspace(0,tmax,n+1)';
    TWanimal.heightsMapped(isnan(TWanimal.heightsMapped)) = [];
    yTW = [0;TWanimal.heightsMapped];
    tSamp = linspace(0,tmax,floor((n+1)/m))';    
else
    m = 1; %#ok<UNRCH>
    n = sum(~isnan(y2));
    tTW = linspace(0,tmax,n+1)';
    y2(isnan(y2)) = [];
    yTW = [0;y2];
    tSamp = linspace(0,tmax,floor((n+1)/m))';
    uDes = [0;u2];
end

nSamps = length(tSamp);
uSamp = zeros(nSamps,2);

x01 = 0;
x02 = 0;
k = 0;
lb = -eps;
ub = 3;

h2 = figure;

for i=2:nSamps
    
    %maximize the force
    x0 = uSamp(i-1,1);
    imax = i;
    tLoop = tTW(1:imax);
    yDesTW = yTW(1:imax);
    
    g1 = @(x)loadMatchingFunc(x,i,tau,a,b,c,d,tSamp,uSamp(:,1),yDesTW,x01);
    
    if g1(x0) == 0
        %do nothing
    else
        x0 = fzero(g1,[lb,ub]);
    end
    
    [uSamp(i,1),~,~,converged] = cornerSearch(g1,x0,p0,lb,ub,100,xTol,fTol,false);
    if ~converged
        error('cornerSearch did not converge. It likely reached the maximum number of iterations. Consider a larger value for p0.')
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    x0 = uSamp(i-1,2);

    g2 = @(x)loadMatchingFunc(x,i,tau,a,b,c,d,tSamp,uSamp(:,2),yDesTW,x02);
    
    if g2(x0) == 0
        %do nothing
    else
        x0 = fzero(g2,[lb,ub]);
    end
    
    [uSamp(i,2),~,~,converged] = cornerSearch(g2,x0,p0,lb,ub,100,xTol,fTol,true);
    if ~converged
        error('cornerSearch did not converge. It likely reached the maximum number of iterations. Consider a larger value for p0.')
    end
    
    figure(h2)
    clf
    subplot(1,2,1)
    plot(tTW(1:i),yDesTW(1:i),'linewidth',2)
    hold on
    [yLo,xLo] = simulateMinusLowpassPL(uSamp(1:i,1),tau,tSamp(1:i),a,b,c,d,0);
    plot(tLoop,yLo,':','linewidth',2)
    [yHi,xHi] = simulateMinusLowpassPL(uSamp(1:i,2),tau,tSamp(1:i),a,b,c,d,0);
    
    x01 = xLo(end);
    x02 = xHi(end);
    
    subplot(1,2,2)
    plot(t,u2,'linewidth',2)
    hold on
    plot(tSamp(1:i),uSamp(1:i,1),'--','linewidth',2)
    plot(tSamp(1:i),uSamp(1:i,2),':','linewidth',2)
    drawnow
end

uInterped = @(u) pchip(tSamp,u,tTW);

lineColors = lines(3);


subplot(1,2,2)
% plot(t,u2,'linewidth',2)
hold on
% plot(tSamp,uSamp(:,1),'--','linewidth',2)
% plot(tSamp,uSamp(:,2),':','linewidth',2)
legend('actual force','lower bound force','upper bound force','location','south','FontSize',8)
title('Applied force','FontSize',8)
xlabel('time (s)','FontSize',8)
ylabel('mN','FontSize',8)
ylim([0,1.5])
box off

subplot(1,2,1)
% plot(tTW,yDesTW,'linewidth',2)
hold on
% plot(tTW,simulateMinusLowpassPL(uInterped(uSamp(:,1)),tau,tTW,a,b,c,d,0),':','linewidth',2)
ylim([0,200])
ylabel('Hz','FontSize',8)
xlabel('time (s)','FontSize',8)
title('Sensory discharge','FontSize',8)
box off

if useExactAnimalRecordings
    sgtitle('(b) Estimated applied force from discharge pattern (animal)','Fontsize',10)
else
    sgtitle('(b) Estimated applied force from discharge pattern (simulation)','Fontsize',10) %#ok<UNRCH>
end

h2.Position(3) = 700;
h2.Position(4) = 200;
set(h2,'renderer','Painters')