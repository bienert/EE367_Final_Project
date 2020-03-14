%Nicole Bienert
%Purpose: Plot MacGregor data and linearly interpolate to grab
%concentration data pts

%inputs:
    %depthSamples: depth samples for the center of each grid in an array
    %numColumns: number of columns
%outputs: Concentrations in units of mol/L

function [Hconcentration,saltConcentration] = mapConcentrations_ver2(depthSamples,numColumns)


%% plot H+  
load SipleDomeHConcentrationMacGregor2007
condY=SipleDomeHConcentrationMacGregor2007(:,2)+5;
condY(length(condY)+1)=3000;
condX=SipleDomeHConcentrationMacGregor2007(:,1);
condX(length(condX)+1)=SipleDomeHConcentrationMacGregor2007(end,1);


%interpolate
vq = interp1(condY,condX,depthSamples);

%plot
figure()
plot(condX,condY)
hold on
plot(vq,depthSamples,'--')

%plot it like macgregor by interleaving to make steps
myLength = length(SipleDomeHConcentrationMacGregor2007(:,1));
condY(1:2:myLength*2)=SipleDomeHConcentrationMacGregor2007(:,2);
condY(2:2:myLength*2-2)=SipleDomeHConcentrationMacGregor2007(2:myLength,2);
condY(myLength*2)=SipleDomeHConcentrationMacGregor2007(end,2);
condY(myLength*2+1)=3000;

condX(2:2:myLength*2)=SipleDomeHConcentrationMacGregor2007(:,1);
condX(1:2:myLength*2)=SipleDomeHConcentrationMacGregor2007(:,1);
condX(myLength*2+1)=SipleDomeHConcentrationMacGregor2007(end,1);
% figure()
plot(condX,condY)

hLegend = legend('Midpts','Linear Interpolation of Midpts','MACGREGOR 2007');
hXlabel = xlabel('Concentration (\mu M)');
hYlabel = ylabel('Depth (m)');
hTitle = title('Siple Dome H^+ Concentration MacGregor 2007');
Aesthetics_Script;


Hconcentration=repmat(vq',1,numColumns)*10^-6;

%% Plot Salt

load SipleDomeSaltConcentrationMacGregor2007
condY=SipleDomeSaltConcentrationMacGregor2007(:,2)+5;
condY(length(condY)+1)=3000;
condX=SipleDomeSaltConcentrationMacGregor2007(:,1);
condX(length(condX)+1)=SipleDomeSaltConcentrationMacGregor2007(end,1);

vq = interp1(condY,condX,depthSamples);
figure()
plot(condX,condY);
hold on
plot(vq,depthSamples,'--');


myLength = length(SipleDomeSaltConcentrationMacGregor2007(:,1));
condY(1:2:myLength*2)=SipleDomeSaltConcentrationMacGregor2007(:,2);
condY(2:2:myLength*2-2)=SipleDomeSaltConcentrationMacGregor2007(2:myLength,2);
condY(myLength*2)=SipleDomeSaltConcentrationMacGregor2007(end,2);
condY(myLength*2+1)=3000;

condX(2:2:myLength*2)=SipleDomeSaltConcentrationMacGregor2007(:,1);
condX(1:2:myLength*2)=SipleDomeSaltConcentrationMacGregor2007(:,1);
condX(myLength*2+1)=SipleDomeSaltConcentrationMacGregor2007(end,1);

% figure()
plot(condX,condY)

hLegend = legend('Midpts','Linear Interpolation of Midpts','MACGREGOR 2007');
hXlabel = xlabel('Concentration (\mu M)');
hYlabel = ylabel('Depth (m)');
hTitle = title('Siple Dome ss Cl^- Concentration MacGregor 2007');
Aesthetics_Script

saltConcentration=repmat(vq',1,numColumns)*10^-6;
