lear
close all
clc

ctr = importdata('CTrTorqueWaveform.csv');

RAHanimal = load('rampAndHold.mat');
RAHanimal.heightsMapped(isnan(RAHanimal.heightsMapped)) = [];
RAHanimal.heightsMapped = [0;RAHanimal.heightsMapped];
n = sum(~isnan(RAHanimal.heightsMapped));

t = ctr.data(:,1);

numSteps = length(t);
tmax = max(t);

t = linspace(0,tmax,n)';
t = t(8:end-9);
y = RAHanimal.heightsMapped(8:end-9);

figure
plot(t,y)

figure
plot(t,y,'+')
ax = gca;
ax.XScale = 'log';
ax.YScale = 'log';

b = log10(y);
A = [log10(t),ones(size(t))];

x = linsolve(A,b);

s = x(1);
a = 10.^x(2);

hold on
plot(t,a*t.^s)

ylim([20,200])
xlim([.1,.7])

b = 1 - 1/s