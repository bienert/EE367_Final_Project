%Nicole Bienert
%Purpose: Compute average conductivity using a CMP and presented in
%winebrenner2004. 

clc; 
close all; 
clear; 

% number of measurements
N = 20*3; %number of measurements per transect
numGridsY=20;
numGridsX=40;

gridResolution=[numGridsY numGridsX];
thickness= 1000;
maxOffset=10000;
hGrid=thickness/numGridsY;
wGrid=maxOffset/numGridsX;

er_ice=3.18; %reletive permittivity of ice
mu=4*pi*10^(-7); %permeability of vacuum
e_0 = 8.854e-12; %permittivity of vaccume


%% temperature options
%Instructions: uncomment whichever distribution you want to test

% temp=fspecial('gaussian', [2*numGridsY numGridsX], 10);
% temp=temp(1:numGridsY,:)./max(max(temp));
% temp=temp*25-30;

% temp=-15*ones(gridResolution);

%no temp anomoly
% temp=repmat(linspace(-26,0,numGridsY)',1,numGridsX);

% %water channel
r = 5; %radius of filter in pixels
[x, y] = meshgrid(-numGridsX/2:numGridsX/2-1,-numGridsY:numGridsY);
b=sqrt(x.^2+y.^2); %distance of each pixel to the cut off
lpf=b<=r; %create low pass filter
temp=lpf(1:numGridsY,:)*15-14;

%% convert temp to cond
[Hconcentration,saltConcentration] = mapConcentrations_ver2_mat([hGrid/2:hGrid:thickness-hGrid/2],numGridsX);
cond=temp2cond_v2_mat(temp,Hconcentration,saltConcentration);

scaling=100;
I=cond*scaling; %cond was too small


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set up start/ end pts for transects
step=(maxOffset/2)/N;

txLoc=fliplr([1:N]*step);
rxLoc=[1:N]*step+maxOffset/2;

% create masks, where \delta r means that the ray passed through that grid
% and 0 is where a ray didn't pass through that grid
for n=1:length(txLoc)
    masks(:,:,n)=pathLenMasks_ver6(txLoc(n),rxLoc(n),wGrid,hGrid,numGridsY,numGridsX,0);
end

%Display masks (ray geometry)
gcf=figure()
for k=1:size(masks,3)
    imagesc(masks(:,:,k))
    title('Ray Geometry')
    pause(0.001)
%     saveas(gcf,['figures\ADMMTV\ADMMTV_p=10_lam=1div10_40x20grids_for_1000x1000m_ThreeOffsetPasses\',num2str(k)],'jpg')
end

N=length(txLoc);


% define function handle for image formation
% please go over them and understand what they do
Afun    = @(x) squeeze(sum(sum(masks .* repmat(x, [1 1 N]),1),2) );

% noise parameter - standard deviation
sigma = 0.0001;

% simulated measurements
atten=-1/2*sqrt(mu./(e_0.*er_ice)).*(Afun(I)+sigma.*randn([N 1]));  %e-az*e-bz=e-(a+b)z

num=atten(2:end)-atten(1);
den=Afun(ones(size(I)));
den=-(den(2:end)-den(1));
line=num./den; 
% figure()
% plot(num)
% figure()
% plot(den)
figure()
plot(line)
title('Normalized attenuation')

%compute average conductivity from attenuation
measuredCond=mean(line(3:end))*2/sqrt(mu./(e_0.*er_ice))/scaling

%calculate PSNR
x=ones(gridResolution)*measuredCond;
I=I/scaling; %rescale for PSNR calculation
MSE = 1/size(I,1)/size(I,2)*sum(sum((x - I).^2));
PSNR = 10*log10(max(I(:).^2)/MSE);
    
%compute temp from conductivity
T=repmat(mean(cond2temp_v3_mat(x,Hconcentration,saltConcentration)),gridResolution(1),1);


figure()
subplot(1,2,1)
imagesc(temp)
colorbar
lim = caxis;
subplot(1,2,2)
imagesc(T)
caxis(lim)
colorbar


figure()
imagesc([0.5 numGridsX-0.5],[0.5 numGridsY-0.5],temp)
colorbar
lim = caxis;
cmocean('thermal')
%plot grids
hold on
% for k = 1:numGridsX-1
%     plot([k k],[0 numGridsY],'Color',[0.6,0.6,0.6])
% end
% hold on
% for k = 1:numGridsY-1
%     plot([0 numGridsX],[k k],'Color',[0.6,0.6,0.6])
% end
hTitle=title('Gaussian Temp Distribution')
hYlabel=ylabel('Depth (m)');
hXlabel=xlabel('Antenna Separation (m)');
%colorbar
%change tick marks
xt=xticks;
xticklabels(wGrid*xt);
yt=yticks;
yticklabels(hGrid*yt);
Aesthetics_Script;


figure()
imagesc([0.5 numGridsX-0.5],[0.5 numGridsY-0.5],T)
colorbar
caxis(lim)
cmocean('thermal')
%plot grids
hold on
% for k = 1:numGridsX-1
%     plot([k k],[0 numGridsY],'Color',[0.6,0.6,0.6])
% end
% hold on
% for k = 1:numGridsY-1
%     plot([0 numGridsX],[k k],'Color',[0.6,0.6,0.6])
% end
hTitle=title({'CMP Inversion';['PSNR=',num2str(PSNR),'dB']})
hYlabel=ylabel('Depth (m)');
hXlabel=xlabel('Antenna Separation (m)');
%colorbar
%change tick marks
xt=xticks;
xticklabels(wGrid*xt);
yt=yticks;
yticklabels(hGrid*yt);
Aesthetics_Script;






