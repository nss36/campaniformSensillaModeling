clear
close all
clc

tau = .03;

tmax = 100;
dt = 1e-3;
tvec = (0:dt:tmax)';
numSteps = length(tvec);

A = 2;
T = 4; %.4;
b = 5;
a = 5;

%Create a ramp-and-hold stimulus with slope A/T and maximum value A.
u = A*tvec/T;
u(u > A) = A;

%Initialize x, the dynamic threshold.
x0 = 0;

%Define x's dynamics and simulate
fPL = @(t,x) 1/tau*sign(interp1(tvec,u,t) - x).*(abs(interp1(tvec,u,t) - x)).^b;

options = odeset('AbsTol',1e-9,'RelTol',1e-6);
[~,x] = ode45( fPL,tvec,x0,options);

%Calculate the model's output, the force (u) relative to a threshold (x)
%scaled by a.
y = a*(u-x);

colorList = lines(3);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOT STEADY-STATE RESPONSE TO RAMP INPUT
% FIGURE S1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

h1 = figure;
plot(tvec,u,'linewidth',2)
hold on
plot(tvec,A/T*tvec - (tau*A/T)^(1/b),'--','linewidth',2,'color',colorList(3,:))
plot(tvec,x,':','linewidth',2,'color',colorList(2,:))
xlim([0,T])
% ylim([0,A+.5])
xlabel('time (s)')
ylabel('activation (n.d.)')
legend('u(t) = A/T\cdott','u(t-\Deltat) = A/T\cdott-(\tau\cdotA/T)^{1/b}','x(t), simulated','location','southeast')
title('x(t) \rightarrow u(t-\Deltat) as time progresses') 
grid on
h1.Position(3) = 400;
h1.Position(4) = 325;
drawnow

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEMONSTRATE BASIC DYNAMICS OF THE MODEL
% FIGURE S2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

h2 = figure;
plot(tvec,u,'linewidth',2)
hold on
plot(tvec,x,':','linewidth',2,'color',colorList(2,:))
plot(tvec,y,'--','linewidth',2,'color',colorList(3,:))
xlim([0,10*T])
% ylim([0,A+.5])
xlabel('time (s)')
ylabel('activation (n.d.)')
legend('u(t)','x(t)',['y(t) = ',num2str(a),'\cdot(u(t)-x(t))'],'location','east')
title({'x follows u, y responds to changes in u,';'and y adapts to tonic u'})
title('') 
grid on
h2.Position(3) = 400;
h2.Position(4) = 325;
drawnow

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% COMPARING ANALYTICAL AND NUMERICAL RESPONSE TO STEP INPUT
% FIGURE S3
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Numerically simulate power-law decay from x0 = A.
x0 = A;
fPL = @(t,x) 1/tau*sign(0 - x).*(abs(0 - x)).^b;
options = odeset('AbsTol',1e-9,'RelTol',1e-6);
[t3,y] = ode45( fPL,tvec,x0,options);

%Analytically calculate power-law decay from x0 = A. These formulas are
%derived in the supplementary materials.
beta = 1-b;
s = 1/beta;
B = (-beta)^(1/beta);
delta_t = tau*(A^(1-b))/(b-1);
f = B*((tvec+delta_t)/tau).^s;

%Plot figure S3
h = figure;
subplot(1,3,1)
plot(tvec,y,'linewidth',2)
hold on
plot(tvec,f,':','linewidth',2)
xlabel('time (s)')
ylabel('activation (n.d.)')
title({'(a) Comparison btwn analytical','and simulation solution'})
legend('simulation','analytical','location','northeast')

subplot(1,3,2)
plot(tvec,y,'linewidth',2)
hold on
plot(tvec,f,':','linewidth',2)
xlabel('time (s)')
xticks(logspace(-3,2,6))
ylabel('activation (n.d.)')
ax = gca;
ax.XScale = 'log';
ax.YScale = 'log';
title({'(b) Comparison with','logarithmic axes'})

subplot(1,3,3)
plot(tvec,abs(f-y),'linewidth',2)
xlabel('time (s)')
ylabel('activation (n.d.)')
xticks(logspace(-3,2,6))
ax = gca;
ax.XScale = 'log';
ax.YScale = 'log';
title({'(c) Error btwn analytical and','simulation solutions'})

h.Position(3) = 700;
h.Position(4) = 250;