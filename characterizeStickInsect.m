% characterizeStickInsect
% Nicholas Szczecinski
% Department of Mechanical and Aerospace Engineering
% West Virginia University
% 2 April 2021

clear
close all
clc

load('stickInsectWalkingForces')
load('stick6AParams')
load('stick6BParams')

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RESPONSE TO SCALED DYNAMIC FORCES
% FIGURE 7
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Estimate how the model responds to scaling of dynamic forces. 
%Compare to Figure 9 of Zill et al. 2021, J Neurophys
peakForceVals = [.57,1.1,1.8,2.3,2.8];
trial = 5;
tmax = .77; %.77; %Time at which Zill et al. 2021 cuts off recording
tmin = .013; %Reduce transient startup signals
maxFfig = figure;
indexShift = 1; %To correlate dudt and y, we may need to shift y by 1 index. This is because it takes time for the system to respond to the rate of change.

tempColors = lines(10);
customColor = [tempColors(2,:);tempColors(1,:);tempColors(7,:);tempColors(3,:);[.5,.5,.5]];

tau6B = stick6BParams.tau;
a6B = stick6BParams.a;
b6B = stick6BParams.b;
c6B = stick6BParams.c;
d6B = stick6BParams.d;

hold on
for j = 1:length(peakForceVals)
    uLoop = -peakForceVals(j)*cleanedData{trial}.u/max(-cleanedData{trial}.u);
    t = cleanedData{trial}.t;
    
    dudt = centeredDiff(t,uLoop);

    u0 = uLoop(1);
    y0 = cleanedData{trial}.sens6B(1);
    y6B = simulateMinusLowpassPL(uLoop,tau6B,t,a6B,b6B,c6B,d6B,u0 - (y0 - c6B*u0 - d6B)/(a6B*1));

    subplot(2,2,1)
    hold on
    plot(t,uLoop,'linewidth',1,'color',customColor(j,:))
    ylabel('Force (mN)')

    subplot(2,2,2)
    hold on
    plot(t,dudt,'linewidth',1,'color',customColor(j,:))
    ylabel('Rate of force (mN/s)')


    subplot(2,2,3)
    hold on
    plot(t,y6B,'linewidth',1,'color',customColor(j,:))
    ylabel('Model discharge (AP/s)')
    ylim([0,300])
    
    subplot(2,2,4)
    hold on
    scatter3(t(t<tmax & t>tmin),circshift(dudt(t<tmax & t>tmin),indexShift),y6B(t<tmax & t>tmin),15,customColor(j,:),'filled')
%     xlim([0,.55])
    ylim([-30,20])
    zlim([0,300])
    zticks(0:50:300)
end



subplot(2,2,1)
title('(a) Single step naturalistic torque waveform','FontSize',8)
xlim([0,tmax])
subplot(2,2,2)
title('(b) Rate of change of waveform','FontSize',8)
xlim([0,tmax])
ylim([-7,15])
subplot(2,2,3)
title('(c) 6B Discharge in response','FontSize',8)
xlim([0,tmax])
subplot(2,2,4)
title('(d) Discharge reflects dF/dt','FontSize',8)
xlabel('t (s)')
ylabel('dF/dt (mN/s)')
zlabel('Model discharge (AP/s)')
view([90 0])
grid on

maxFfig.Position(3) = 700;
