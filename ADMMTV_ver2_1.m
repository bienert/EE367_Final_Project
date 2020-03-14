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
    %ver1_6: initialize matricies to be 0
%ver2: also invert for reflection coeffecient

clc; 
close all; 
clear; 

% number of measurements
N = 20; %number of measurements per transect
numGridsY=20;
numGridsX=40;
p1 = 50; 
p2 = 100; 

lambda1 = 10; 
lambda2=0.08;
numIters    = 200;  % number of ADMM iterations


thickness= 1000;
maxOffset=1000;
hGrid=thickness/numGridsY;
wGrid=maxOffset/numGridsX;

% figure()
% I=fspecial('gaussian', [2*numGridsY numGridsX], 10);
% I=I(1:numGridsY,:);

%no temp anomoly
% I=repmat([1:numGridsY]',1,numGridsX);

% %water channel
r = 5; %radius of filter in pixels
[x, y] = meshgrid(-numGridsX/2:numGridsX/2-1,-numGridsY:numGridsY);
b=sqrt(x.^2+y.^2); %distance of each pixel to the cut off
lpf=b<=r; %create low pass filter
cond=lpf(1:numGridsY,:);

condTrue=lpf(1:numGridsY,:);

reflTrue= zeros(size(cond));
reflTrue(end,:)=1;

I=[condTrue,reflTrue]; %real image we are trying to invert for


gridResolution=[numGridsY numGridsX];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set up start/ end pts for transects
%transect 1
txLoc1=400;
step=(maxOffset-txLoc1)/N;
rxLoc=txLoc1+[1:N]*step;
%transect 2
txLoc2= 700;
step=(txLoc2)/N;
rxLoc=[rxLoc,txLoc2-[1:N]*step];
%transect 3
txLoc3=0;
step=(maxOffset-txLoc3)/N;
rxLoc=[rxLoc,txLoc3+[1:N]*step];
txLoc=[repmat(txLoc1,1,N),repmat(txLoc2,1,N),repmat(txLoc3,1,N)];


% create masks, where \delta r means that the ray passed through that grid
% and 0 is where a ray didn't pass through that grid
for n=1:length(txLoc)
    rayGeometry(:,:,n)=pathLenMasks_ver5(txLoc(n),rxLoc(n),wGrid,hGrid,numGridsY,numGridsX,0);
end


%create reflection coeffecient masks and concatinate with ray geometry
%if reflection hits between two grids, then we round down
for n=1:length(txLoc)
    reflMask(:,:,n)=reflMasksFun_ver1_1(txLoc(n),rxLoc(n),wGrid,hGrid,numGridsY,numGridsX,0);
end

% %Display masks (ray geometry)
% gcf=figure()
% for k=1:size(rayGeometry,3)
%     imagesc(rayGeometry(:,:,k)+max(max(rayGeometry(:,:,k)))*2*reflMask(:,:,k))
%     title('Geometry')
% %     saveas(gcf,['figures\ADMMTV\ADMMTV_p=10_lam=1div10_40x20grids_for_1000x1000m_ThreeOffsetPasses\',num2str(k)],'jpg')
% 
% %     saveas(gcf,['figures\ADMMTV\ADMMTV_p=10_lam=1div10_40x20grids_for_1000x1000m_ThreeOffsetPasses\ADMMTV_p=',num2str(p),'_lam=1div',num2str(1/lambda),'_',num2str(numGridsX),'x',num2str(numGridsY),'grids_for_',num2str(maxOffset),'x',num2str(thickness),'m_a_ray_geometry_',num2str(k)],'jpg')
%     pause(0.1)
% end

N=length(txLoc);

M1=[ones(gridResolution),zeros(gridResolution)]; %mask to isolate conductivity
M2=[zeros(gridResolution),ones(gridResolution)]; %mask to isolate reflection coeffecient

G=[rayGeometry,reflMask];
imageResolution=[size(G,1) size(G,2)];
% define function handle for image formation
% please go over them and understand what they do
Afun    = @(x) squeeze(sum(sum(G .* repmat(x, [1 1 N]),1),2) );
Atfun   = @(x) sum(repmat(reshape(x, [1 1 N]), [imageResolution(1) imageResolution(2) 1]).*G, 3);
opDtDxfun = @(x) opDtx(opDx(reshape(x,[imageResolution(1) imageResolution(2)])));
Atilfun  = @(x) reshape(Atfun(Afun(reshape(x,[imageResolution(1) imageResolution(2)])))+p1.*M1.*opDtDxfun(M1.*reshape(x,[imageResolution(1) imageResolution(2)])) +p2.*M2.*reshape(x,[imageResolution(1) imageResolution(2)]),[prod(imageResolution) 1]);

% noise parameter - standard deviation
sigma = 0.01;

% simulated measurements
b   = Afun(I) + sigma.*randn([N 1]);

%% %%%%%%%%%%%% ADMM%%%%%%%%%%%%%%
% ADMM with (anisotropic) TV
rowsCols = [1 2];

x = randn(imageResolution);
z1 = zeros([imageResolution(1) imageResolution(2) 2]);
u1 = zeros([imageResolution(1) imageResolution(2) 2]);
z2 = zeros([imageResolution(1) imageResolution(2) 2]);
u2 = zeros([imageResolution(1) imageResolution(2) 2]);

M1=[ones(gridResolution),zeros(gridResolution)]; %mask to isolate conductivity
M2=[zeros(gridResolution),ones(gridResolution)]; %mask to isolate reflection coeffecient

gcf=figure()
subplot(1,2,1)
imagesc(I)
subplot(1,2,2)
imagesc(x)
pause(0.5)
%     saveas(gcf,['figures\ADMMTV\ADMMTV_p=10_lam=1div10_40x20grids_for_1000x1000m_ThreeOffsetPasses\',num2str(k+N*3+1)],'jpg')

%     saveas(gcf,['figures\ADMMTV\ADMMTV_p=10_lam=1div10_40x20grids_for_1000x1000m_ThreeOffsetPasses\ADMMTV_p=',num2str(p),'_lam=1div',num2str(1/lambda),'_',num2str(numGridsX),'x',num2str(numGridsY),'grids_for_',num2str(maxOffset),'x',num2str(thickness),'m_b_inversion_',num2str(0)],'jpg')

for ind=1:numIters

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % x update
    v=z1-u1;
    q=z2-u2;
    btil=reshape(Atfun(b)+p1*M1.*opDtx(v)+p2*M2.*(q(:,:,1)+q(:,:,2)),[prod(imageResolution) 1]);%b is a matrix of size x, but for pcg it needs to be a column vector
    x=pcg(Atilfun,btil,1e-12,25,[],[],x(:));
    x=reshape(x,[imageResolution(1) imageResolution(2)]);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % z1 update - soft shrinkage

    kappa1 = lambda1/p1;
    v = (opDx(M1.*x)+u1); 
    z1 = max(1 - kappa1/abs(v),0) .* v;
    
        % z2 update - soft shrinkage

    kappa2 = lambda2/p2;
    w = ((M2.*x)+u2); 
    z1 = max(1 - kappa2/abs(w),0) .* w;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % u1 update    
    u1 = u1 + opDx(M1.*x)-z1;
    
    % u2 update    
    u2 = u2 + M2.*x-z2;


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % PSNR and residuals
    MSE = 1/size(I,1)/size(I,2)*sum(sum((x - I).^2));
    PSNR(ind) = 10*log10(max(I(:).^2)/MSE);
    
    term1=b-Afun(x);
    term2=opDx(x);
    residual(ind)=0.5*sum(term1(:).^2) + lambda1.*sum( abs(term2(:)) );  
    
    
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



