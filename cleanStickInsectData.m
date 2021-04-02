% tuneStickInsect
% Nicholas Szczecinski
% Department of Mechanical and Aerospace Engineering
% West Virginia University
% 2 April 2021

%This script does the following:
%-Because all force stimuli are saved as 'positive' forces but labeled as
%flexion or extension, we must convert them into positive and negative
%forces for the model.
%-To aid the model tuning, delay the sensory recording by one 'time step'.
%-Because large and small sensilla APs could be identified in 6B recordings
%but not from 6A recordings, we need to separately identify them from the 

close all
clear
clc

trials = {  'trial1.csv';...
            'trial2.csv';...
            'trial3.csv';...
            'trial4.csv';...
            'trial5.csv';...
            'trial6.csv';...
            'trial7.csv';...
            'trial8.csv'};
        

%Plot all of the forces together, as a sanity check.

hForce = figure;
title('walking force stimuli')
hold on

trial = cell(8,1);

peakForceVals = [.57,1.1,1.8,2.3,2.8];

for i=1:8
    %Import a trial and extract some info from it.
    trial{i} = importdata(trials{i});
    t = trial{i}.data(:,1)/100; %percent stance
    dt = mean(diff(t));
    t = [0;t+dt/2];
    u = trial{i}.data(:,4); %force
    dudt = [0;trial{i}.data(:,5)]; %df/dt
    sens6A = [0;trial{i}.data(:,3)];
    
    %These are flexor forces, so the sign of the force must be flipped, and
    %we need to retrieve column 6 from the raw data for 6B discharge (not
    %column 2 as for the extensor forces).
    if i >= 5
        %here, the flexor torque is recorded, so invert the sign of the force.
        u = -u;
        sens6B = [0;trial{i}.data(:,6)]; %6B, all sensilla
        sens6Blarge = [0;trial{i}.data(:,2)]; %6B, large sensilla only
    else
        sens6B = [0;trial{i}.data(:,2)]; %6B, all sensilla, but probably only large sensilla due to difficulty of recording.
    end
    u(end+1) = 0; %#ok<SAGROW>
    
    %Now that we've cleaned up this data to suit our needs, let's overwrite
    %the old data.
    cleanedData{i}.t = t; %#ok<SAGROW>
    cleanedData{i}.sens6A = sens6A; %#ok<SAGROW>
    cleanedData{i}.sens6B = sens6B; %#ok<SAGROW>
    cleanedData{i}.u = u; %#ok<SAGROW>
    if i >= 5
        %For the flexor torques, we have access to the discharge of large
        %6B caps, only.
        cleanedData{i}.sens6Blarge = sens6Blarge; %#ok<SAGROW>
    end
    
    plot(t,u,'linewidth',1)
end

figure(hForce)
xlabel('time (s)')
ylabel('force (mN)')

save('stickInsectWalkingForces.mat','cleanedData');
