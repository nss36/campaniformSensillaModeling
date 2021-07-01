% Effect of cocontraction on CS load encoding
% Nicholas Szczecinski
% West Virginia University
% 25 May 2021

rExt = 1; %mm
rFlx = 1; %mm
Ltibia = 20; %mm

%Position of the CS on the dorsal surface of the tibia. x = 0 is the FTi
%joint. Increasing values of x are distal.
rCS = 1.1; %mm

%Set the number of cocontractions to test. Establish a vector c of
%fractions from 0 to 1, excluding 1.
numCoContract = 10;
c = linspace(0,(numCoContract-1)/numCoContract,numCoContract);
% c = linspace((numCoContract-1)/numCoContract,0,numCoContract);

%Set the number of loads applied perpendicular to the tibia at the TiTar
%joint to test.
numLoad = 5;
maxLoad = 4; %Positive is directed inward, negative is directed outward.
Fload = linspace(0,maxLoad,numLoad); %mN

%For each load and degree of cocontraction, we will calculate the extensor
%force, the flexor force, and the force the condyle applies to keep the FTi
%joint stationary. Each index corresponds to the conditions (load and
%cocontraction).
Fext = NaN(numLoad,numCoContract); %mN
Fflx = NaN(numLoad,numCoContract); %mN
FO = NaN(numLoad,numCoContract); %mN

loadInds = 3;
cocontractInds = 1:numCoContract/2;

colors = zeros(length(cocontractInds),3);
colors(:,1) = linspace(0,1,length(cocontractInds));
beamDiagrams = figure;
lw = 1;

%Simultaneously solve the sum of forces in the y direction, moments about
%the FTi joint, and the relative cocontraction of the flexor and extensor.
for i=1:numLoad
    for j=1:numCoContract
        A = [1,1,-1;rExt,-rFlx,0;c(j),-1,0];
        b = [-Fload(i);Ltibia*Fload(i);0];
        x = linsolve(A,b);
        
        Fext(i,j) = x(1);
        Fflx(i,j) = x(2);
        FO(i,j) = x(3);
        
        if any(loadInds == i) && any(cocontractInds == j)
            figure(beamDiagrams);
            hold on
            subplot(3,1,1)
            plot(-rExt*[1,1],[0,-Fext(i,j)],'color',colors(j,:),'linewidth',lw)
            hold on
            plot(0*[1,1],[0,FO(i,j)],'color',colors(j,:),'linewidth',lw)
            plot(rFlx*[1,1],[0,-Fflx(i,j)],'color',colors(j,:),'linewidth',lw)
            plot(Ltibia*[1,1],[0,-Fload(i)],'color',colors(j,:),'linewidth',lw)
            grid on
            
            subplot(3,1,2)
            hold on
            xVals = [-rExt,-rExt,0,0,rFlx,rFlx,Ltibia,Ltibia];
            V = [0,-Fext(i,j),-Fext(i,j),FO(i,j)-Fext(i,j),FO(i,j)-Fext(i,j),FO(i,j)-Fext(i,j)-Fflx(i,j),FO(i,j)-Fext(i,j)-Fflx(i,j),0];
            plot(xVals,V,'color',colors(j,:),'linewidth',lw)
            grid on
            
            subplot(3,1,3)
            hold on
            xVals = [-rExt,0,rFlx,Ltibia];
            M = [0,-rExt*Fext(i,j),-(Ltibia - rFlx)*Fload(i),0];
            plot(xVals,M,'color',colors(j,:),'linewidth',lw)
            grid on
        end
    end
end

figure(beamDiagrams)
subplot(3,1,1)
hold on
% legend(
    
if rCS < 0
    Mcs = (-rExt - rCS)*Fext;
elseif rCS < rFlx
%     Mcs = Fext*(1 - rCS/(rExt + rFlx));
    Mcs = rCS/rFlx*(Ltibia - rFlx)*(FO - Fext - Fflx) - (1-rCS/rFlx)*rExt*Fext;
else
    Mcs = (Ltibia - rCS)*(FO - Fext - Fflx);
end

[FLOAD,C] = meshgrid(Fload,c);
FLOADt = FLOAD';
Ct = C';

figure
scatter3(FLOADt(:),Ct(:),-Mcs(:))
xlabel('Load at tarsus, perpendicular to tibia (mN)')
ylabel('Degree of cocontraction')
zlabel('Bending moment at CS location (\muNm)')