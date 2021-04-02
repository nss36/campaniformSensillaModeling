% tuneStickInsect
% Nicholas Szczecinski
% Department of Mechanical and Aerospace Engineering
% West Virginia University
% 2 April 2021

close all
clear
clc

%Load "cleanedData" struct, which includes t, u, 6A, and 6B.
load('stickInsectWalkingForces.mat');

%Pick one extensor torque and one flexor torque trial to use to tune the
%parameter values for the 6A and 6B models, respectively.
extTrialToTune = 1;
flxTrialToTune = 5;

%make u for tuning 6A
t = cleanedData{extTrialToTune}.t;
u6A = cleanedData{extTrialToTune}.u;
sens6A = cleanedData{extTrialToTune}.sens6A;

%make u for tuning 6B
t = cleanedData{flxTrialToTune}.t;
u6B = -cleanedData{flxTrialToTune}.u;
sens6B = cleanedData{flxTrialToTune}.sens6B;

%make u for tuning 6B large. We need to use one of the ext trials, because
%these recordings only have large 6B action potentials.
u6Blarge = -cleanedData{flxTrialToTune}.u; %force
sens6Blarge = cleanedData{flxTrialToTune}.sens6Blarge; %6B large

%Identify the initial stimulus value, which will be used to initialize x,
%the dynamic threshold.
u06A = u6A(1);
u06B = u6B(1);
u06Blarge = u6Blarge(1);
y06A = sens6A(1);
y06B = sens6B(1);
y06Blarge = sens6Blarge(1);

%Define functions that simulate the response of each CS group, and return
%it as a vector with the same length as sens6A, sens6B, and sens6Blarge.
fsim6A = @(x,u) simulateMinusLowpassPL(u,x(1)/1000,t,x(2)*1000,x(3),x(4),x(5),u06A - (y06A - x(4)*u06A - x(5))/(x(2)*1000));
fsim6B = @(x,u) simulateMinusLowpassPL(u,x(1)/1000,t,x(2)*1000,x(3),x(4),x(5),u06B - (y06B - x(4)*u06B - x(5))/(x(2)*1000));
fsim6Blarge = @(x,u) simulateMinusLowpassPL(u,x(1)/1000,t,x(2)*1000,x(3),x(4),x(5),u06Blarge - (y06Blarge - x(4)*u06Blarge - x(5))/(x(2)*1000));

%Using fsim functions, define cost functions that calculate the mean
%squared error between each recording and the corresponding model output.
f6A = @(x)sum( (fsim6A(x,u6A) - sens6A).^2 );
f6B = @(x)sum( (fsim6B(x,u6B) - sens6B).^2 );
f6Blarge = @(x)sum( (fsim6Blarge(x,u6Blarge) - sens6Blarge).^2 );

% %Approximate initial values based on Zill, Bueschges, and Schmitz 2011.
% tauA = 3.543; %This is in ms. Note that the "fsim" functions below scale this value by multiplying by 1e-3.
% aA = 1.024;
% bA = 2.185;
% cA = 12; %From Zill, Buescghes, and Schmitz 2011, Fig. 4h
% dA = -9; %From Zill, Buescghes, and Schmitz 2011, Fig. 4h
% 
% tauB = 3.8; %3.543; %This is in ms. Note that the "fsim" functions below scale this value by multiplying by 1e-3.
% aB = 0.7; %600; 
% bB = 2.3; %1; 
% cB = 54; %22; %From Zill, Buescghes, and Schmitz 2011, Fig. 4g
% dB = -41; %15; %From Zill, Buescghes, and Schmitz 2011, Fig. 4g

%Initial guesses for 6A parameter values.
% x06A = [tau;a/3;b;cA;dA];
hResults6A = figure;
subplot(4,1,1)
plot(t,u6A,'linewidth',1)
subplot(4,1,2)
plot(t,sens6A,'linewidth',1)
subplot(4,1,3)
drawnow

% x06B = [tau;2*a/3;b;cB;dB];
hResults6B = figure;
subplot(4,1,1)
plot(t,u6B,'linewidth',1)
subplot(4,1,2)
plot(t,sens6B,'linewidth',1)
subplot(4,1,3)
ylabel('initial')
drawnow

% x06Blarge = [tau;a*3;1;0;-20];
hResults6Blarge = figure;
subplot(4,1,1)
plot(t,u6Blarge,'linewidth',1)
hold on
plot(t,centeredDiff( t,u6Blarge))
subplot(4,1,2)
plot(t,sens6B,'linewidth',1)
hold on
plot(t,sens6Blarge,'linewidth',1)
subplot(4,1,3)
plot(t,sens6Blarge,':','linewidth',1)
ylabel('initial')
drawnow

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% USE GENETIC ALGORITHM TO FIT RESPONSES.
% RUNS RAPIDLY AND GIVES GOOD RESULTS.

lb = [.1;0;1;0;-100];
lbLarge = lb;
ub = [10;2;7;100;100];

maxGen = 1000;
convCrit = maxGen;
popSize = 200;
crossFrac = []; %0.5
crossMeth = 'single';
mutProb = [];
seed = 1;
toPar = true;
toPlot = 1;
toPrint = 2;

x06A = GA( f6A,[lb,ub],[],convCrit,popSize,crossFrac,crossMeth,mutProb,seed,toPar,toPlot,toPrint);
x06B = GA( f6B,[lb,ub],[],convCrit,popSize,crossFrac,crossMeth,mutProb,seed,toPar,toPlot,toPrint);
x06Blarge = GA( f6B,[lbLarge,ub],[],convCrit,popSize,crossFrac,crossMeth,mutProb,seed,toPar,toPlot,toPrint);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% USE GRADIENT-BASED SEARCH TO FIT RESPONSES.
% RUNS RAPIDLY AND GIVES GOOD RESULTS.
optionsfmincon = optimoptions('fmincon','OptimalityTolerance',1e-12,'StepTolerance',1e-12,'UseParallel',false,'Display','iter-detailed');
xf6A = fmincon( f6A,x06A,[],[],[],[],lb,ub,[],optionsfmincon);
xf6B = fmincon( f6B,x06B,[],[],[],[],lb,ub,[],optionsfmincon);
xf6Blarge = fmincon( f6Blarge,x06Blarge,[],[],[],[],lbLarge,ub,[],optionsfmincon);
    
% Values used in manuscript:
%     xf6A = [
% 
%         9.6782;
%         0.2650;
%         1.6750;
%        17.7500;
%       -22.5000];
% 
%     xf6B = [
% 
%         1.6591;
%         0.6050;
%         3.3251;
%         5.7500;
%        10.4999];
% 
%     xf6Blarge = [
% 
%         1.6488;
%         0.6132;
%         3.3103;
%         5.6594;
%       -99.4475];

figure(hResults6A)
hold on
subplot(4,1,4)
plot(t,fsim6A(xf6A,u6A),'linewidth',1)
hold on
plot(t,sens6A,':','linewidth',1)
ylabel('final')

figure(hResults6B)
hold on
subplot(4,1,4)
plot(t,fsim6B(xf6B,u6B),'linewidth',1) %#ok<*UNRCH>
hold on
plot(t,sens6B,':','linewidth',1)
ylabel('final')

figure(hResults6Blarge)
hold on
subplot(4,1,4)
plot(t,fsim6Blarge(xf6Blarge,u6Blarge),'linewidth',1) %#ok<*UNRCH>
hold on
plot(t,sens6Blarge,':','linewidth',1)
ylabel('final')

stick6AParams.tau = 1e-3*xf6A(1);
stick6AParams.a = 1e3*xf6A(2);
stick6AParams.b = xf6A(3);
stick6AParams.c = xf6A(4);
stick6AParams.d = xf6A(5);

stick6BParams.tau = 1e-3*xf6B(1);
stick6BParams.a = 1e3*xf6B(2);
stick6BParams.b = xf6B(3);
stick6BParams.c = xf6B(4);
stick6BParams.d = xf6B(5);

save('stick6AParams.mat','stick6AParams');
save('stick6BParams.mat','stick6BParams');

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOT TUNED MODEL RESPONSES
% PLOTS FROM FIGURE 6
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


filename = 'comparison6A6B.xls';

columns = {'A','H','O','V','AC','AJ','AQ','AX'};
for i=1:8
    fsim6A = @(x,u) simulateMinusLowpassPL(u,x(1)/1000,t,x(2)*1000,x(3),x(4),x(5),u06A - (y06A - x(4)*u06A - x(5))/(x(2)*1000));
    fsim6B = @(x,u) simulateMinusLowpassPL(u,x(1)/1000,t,x(2)*1000,x(3),x(4),x(5),u06B - (y06B - x(4)*u06B - x(5))/(x(2)*1000));
    
    sim6A = fsim6A(xf6A,cleanedData{i}.u); %%GOTTA MAKE NEW fsim WITH NEW ICS
    sim6B = fsim6B(xf6B,-cleanedData{i}.u); %%GOTTA MAKE NEW fsim WITH NEW ICS
    
    hComparison = figure;
    subplot(3,1,1)
    if i >= 5
        u = -cleanedData{i}.u;
    else
        u = cleanedData{i}.u;
    end
    
    plot(cleanedData{i}.t,u,'linewidth',1)
    ylabel('Stimulus')
    
    subplot(3,1,2)
    plot(cleanedData{i}.t,cleanedData{i}.sens6A,'linewidth',1)
    hold on
    plot(cleanedData{i}.t,sim6A,':','linewidth',2)
    
    ylabel('6A discharge')
    legend('animal','model')
    ylim([0,300])
    
    subplot(3,1,3)
    plot(cleanedData{i}.t,cleanedData{i}.sens6B,'linewidth',1)
    hold on
    plot(cleanedData{i}.t,sim6B,':','linewidth',2)
    
    ylabel('6B discharge')
    legend('animal','model'); %,'model, L')
    xlabel('time (s)')
    ylim([0,300])
    
    hdrs = {'time','force','sens6A','sim6A','sens6B','sim6B'};
    datCell = num2cell([cleanedData{i}.t,cleanedData{i}.u,cleanedData{i}.sens6A,sim6A,cleanedData{i}.sens6B,sim6B]);
    
    writecell([hdrs;datCell],filename,'Range',[columns{i},'1'])
end