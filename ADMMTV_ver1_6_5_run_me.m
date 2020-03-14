%Nicole Bienert
%2/28/2020
%Purpose: Use ADMM with a TV prior to invert for conductivity 

%version history:
%ver1: only invert for conductivity and don't use real numbers
    %ver1_2: uses two transects, one moving to the left and the other to the right
    %      plots v and z in intermediate steps - testing on stanford logo 
    %ver1_3: uses two transects
    %ver1_4: use two transects and select their location
    %ver1_5: use three transects and select their locations
    %ver1_6: initialize matricies to be 0. The old ray mask was not
    %accounting for rays that pass through multiple grids.
        %ver1_6_4: Takes temp as input and inverts for temp. This is the
        %complted version
        %ver1_6_5: uses a 10km transect. Include options for different numbers of
        %transects
%ver2: also invert for reflection coeffecient

clc; 
close all; 
clear; 

N = 20; %number of measurements per transect
numGridsY=20; %number of grids per column
numGridsX=40; %number of grids per row
p=100; %rho controls ADMM step size
lambda=10; %lambda controls the TV prior influence
numIters    = 400;  % number of ADMM iterations
thickness= 1000; %ice sheet thickness in meters
maxOffset=10000; %largest antenna separation in meters

%vars from eqs
hGrid=thickness/numGridsY; %height in meters of a grid
wGrid=maxOffset/numGridsX; %width in meters of a grid
gridResolution=[numGridsY numGridsX]; %grid resolution

%% temperature options
%Instructions: uncomment whichever distribution you want to test

% %gaussian distribution
% temp=fspecial('gaussian', [2*numGridsY numGridsX], 10);
% temp=temp(1:numGridsY,:)./max(max(temp));
% temp=temp*25-30;

%Constant temperature
% temp=10*ones(gridResolution);

%no temp anomoly: linear gradient
% temp=repmat(linspace(-26,0,numGridsY)',1,numGridsX);

% %water channel
r = 5; %radius of filter in pixels
[x, y] = meshgrid(-numGridsX/2:numGridsX/2-1,-numGridsY:numGridsY);
b=sqrt(x.^2+y.^2); %distance of each pixel to the cut off
lpf=b<=r; %create low pass filter
temp=lpf(1:numGridsY,:)*15-14;

%% Transect options
%Instructions: uncomment whichever option you want for transects. You will
%need to modify p and lambda accordingly
% % one transect
% %for one transect, use:
% %     p=4980;
% %     lambda=205;
% % step=(maxOffset)/N;
% % txLoc=0;
% % rxLoc=[txLoc+[1:N]*step];
% % txLoc=zeros(1,length(rxLoc));

% % two transects
% %for two transects, use p=100 and lambda = either 5 or 10
% step=(maxOffset)/N;
% txLoc1=0;
% rxLoc1=[txLoc1+[1:N]*step];
% txLoc2=10000;
% rxLoc2=txLoc2-[1:N]*step;
% txLoc=[ones(1,length(rxLoc1))*txLoc1,ones(1,length(rxLoc2))*txLoc2];
% rxLoc=[rxLoc1,rxLoc2];


% three asymmetric transects
%For hyperparameters, use:
%     p=100; 
%     lambda=10;
%transect 1
txLoc1=4000;
step=(maxOffset-txLoc1)/N;
rxLoc=txLoc1+[1:N]*step;
%transect 2
txLoc2= 7000;
step=(txLoc2)/N;
rxLoc=[rxLoc,txLoc2-[1:N]*step];
%transect 3
txLoc3=0;
step=(maxOffset-txLoc3)/N;
rxLoc=[rxLoc,txLoc3+[1:N]*step];
txLoc=[repmat(txLoc1,1,N),repmat(txLoc2,1,N),repmat(txLoc3,1,N)];

%% convert temp to conductivity
%sample acid and salt concentrations at center of grids
[Hconcentration,saltConcentration] = mapConcentrations_ver2_mat([hGrid/2:hGrid:thickness-hGrid/2],numGridsX);
%use chemistry and temperature info to compute conductivity
cond=temp2cond_v2_mat(temp,Hconcentration,saltConcentration);

scaling=100;
I=cond*scaling; %cond was too small

N=length(txLoc);

%% Masks
%create masks, where \delta r means that the ray passed through that grid
% and 0 is where a ray didn't pass through that grid
for n=1:length(txLoc)
    masks(:,:,n)=pathLenMasks_ver6(txLoc(n),rxLoc(n),wGrid,hGrid,numGridsY,numGridsX,0);
end

%Display masks (ray geometry)
gcf=figure()
for k=1:size(masks,3)
    imagesc(masks(:,:,k))
    title('Ray Geometry')
%     saveas(gcf,['figures\ADMMTV\ADMMTV_p=10_lam=1div10_40x20grids_for_1000x1000m_ThreeOffsetPasses\',num2str(k)],'jpg')

%     saveas(gcf,['figures\ADMMTV\ADMMTV_p=10_lam=1div10_40x20grids_for_1000x1000m_ThreeOffsetPasses\ADMMTV_p=',num2str(p),'_lam=1div',num2str(1/lambda),'_',num2str(numGridsX),'x',num2str(numGridsY),'grids_for_',num2str(maxOffset),'x',num2str(thickness),'m_a_ray_geometry_',num2str(k)],'jpg')
    pause(0.1)
end
 
% figure()
% imagesc(masks(:,:,20))
% cmocean('thermal')
% hold on
% hTitle=title({'Ray Geometry';'Transect 1';''})
% hYlabel=ylabel('Depth (m)');
% hXlabel=xlabel('Antenna Separation (m)');
% colorbar
% %change tick marks
% xt=xticks;
% xticklabels(wGrid*xt);
% yt=yticks;
% yticklabels(hGrid*yt);
% Aesthetics_Script;



%% Function Definitions
% define function handle for image formation
% please go over them and understand what they do
Afun    = @(x) squeeze(sum(sum(masks .* repmat(x, [1 1 N]),1),2) );
Atfun   = @(x) sum(repmat(reshape(x, [1 1 N]), [gridResolution(1) gridResolution(2) 1]).*masks, 3);
opDtDxfun = @(x) opDtx(opDx(reshape(x,[gridResolution(1) gridResolution(2)])));
Atilfun  = @(x) reshape(Atfun(Afun(reshape(x,[gridResolution(1) gridResolution(2)])))+p.*opDtDxfun(x) ,[prod(gridResolution) 1]);

% noise parameter - standard deviation
sigma = 0.0001;

% simulated measurements
b   = Afun(I) + sigma.*randn([N 1]);

%% %%%%%%%%%%%% ADMM%%%%%%%%%%%%%%
%Initialize variables
x = zeros(gridResolution);
z = zeros([gridResolution(1) gridResolution(2) 2]);
u = zeros([gridResolution(1) gridResolution(2) 2]);

gcf=figure()
subplot(1,2,1)
imagesc(I)
subplot(1,2,2)
imagesc(x)
pause(0.5)
%     saveas(gcf,['figures\ADMMTV\ADMMTV_p=10_lam=1div10_40x20grids_for_1000x1000m_ThreeOffsetPasses\',num2str(k+N*3+1)],'jpg')
for ind=1:numIters

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % x update
    v=z-u;
    btil=reshape(Atfun(b)+p*opDtx(v),[prod(gridResolution) 1]);%b is a matrix of size x, but for pcg it needs to be a column vector
    x=pcg(Atilfun,btil,1e-20,25,[],[],x(:));
    x=reshape(x,[gridResolution(1) gridResolution(2)]);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % z update - soft shrinkage

    kappa = lambda/p;
    v = (opDx(x)+u); 
    z = max(1 - kappa/abs(v),0) .* v;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % u update    
    u = u + opDx(x)-z;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % PSNR and residuals
    MSE = 1/size(I,1)/size(I,2)*sum(sum((x - I).^2));
    PSNR(ind) = 10*log10(max(I(:).^2)/MSE);
    
    term1=b-Afun(x);
    term2=opDx(x);
    residual(ind)=0.5*sum(term1(:).^2) + lambda.*sum( abs(term2(:)) );  
    
    
    % plot
    gcf=figure(1)
    subplot(1,3,1)
    imagesc(I)
    title('Actual Cond')
    subplot(1,3,2)
    imagesc(x)
    title('Inverted Cond')
    subplot(1,3,3)
    plot(residual)
    title('Residual')
%   saveas(gcf,['figures\ADMMTV\ADMMTV_p=10_lam=1div10_40x20grids_for_1000x1000m_ThreeOffsetPasses',num2str(ind+N+1)],'jpg')

%   saveas(gcf,['figures\ADMMTV\ADMMTV_p=10_lam=1div10_40x20grids_for_1000x1000m_ThreeOffsetPasses\ADMMTV_p=',num2str(p),'_lam=1div',num2str(1/lambda),'_',num2str(numGridsX),'x',num2str(numGridsY),'grids_for_',num2str(maxOffset),'x',num2str(thickness),'m_b_inversion_',num2str(ind)],'jpg')
%     pause(0.5)
end


x=x/scaling;
T=cond2temp_v3_mat(x,Hconcentration,saltConcentration);
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
subplot(1,2,1)
plot(residual)
hTitle=title('Residual')
hYlabel=ylabel('0.5|Cx-b|_2^2+\lambda|\nabla x|_1');
hXlabel=xlabel('Iteration');
Aesthetics_Script
subplot(1,2,2)
plot(PSNR)
hTitle=title('PSNR')
hYlabel=ylabel('dB');
hXlabel=xlabel('Iteration');
Aesthetics_Script


figure()
imagesc([0.5 numGridsX-0.5],[0.5 numGridsY-0.5],temp)
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
hTitle=title('Subglacial Channel')
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
hTitle=title({'ADMM Inversion';['PSNR=',num2str(PSNR(ind)),'dB']})
hYlabel=ylabel('Depth (m)');
hXlabel=xlabel('Antenna Separation (m)');
%colorbar
%change tick marks
xt=xticks;
xticklabels(wGrid*xt);
yt=yticks;
yticklabels(hGrid*yt);
Aesthetics_Script;



