%Nicole Bienert
%Date: 2/24/20

%Purpose: To create path length masks
%ver2: Calculates ray intercepts and grids the ray passed through for the
%ray going either up or down.
%ver3: Calculates the total ray path
%ver5: accepts start and stop locations instead of just antenna separation
%      uses uniquetol instead of unique to handle issues with floats
%ver6: improved plotting and accounts for a ray passing through the same
%grid 2x

%Notes: grids are indexed starting at 0 in x, but matlab starts at 1, so it a
%ray is supposed to go through x/wGrid=16, then it goes through column 17
%in the matrix. This is accomplished with the ceil operator

%based on: https://www.mathworks.com/matlabcentral/answers/230155-how-to-determine-which-grid-cells-a-line-segment-passes-through
%x and y are numbered like this:
% glacier surface
% (1,20)  (2,20)  (3,20)  ...  (numGridsX,20)
% (1,19)  (2,19)  (3,19)  ...  (numGridsX,20)
% (1,18)    ...     ...   ...     ...

function myMat = pathLenMasks_ver5(txLoc,rxLoc,wGrid,hGrid,numGridsY,numGridsX,display)

startPt=min(txLoc,rxLoc);
endPt=max(txLoc,rxLoc);
% numGridsY=20;
% numGridsX=10;
% hGrid=1000/numGridsY;
% wGrid=1000/numGridsX;
% 
% offset = 1000;
endGrid=endPt/wGrid; %number of grids for this antenna separation
startGrid=startPt/wGrid; %Transmitter location normalized to grid width
thickness=hGrid*numGridsY;
m=numGridsY/(endGrid-startGrid)*2;           % Slope (or slope array)        
b=numGridsY-m*(endGrid);                       % Intercept (or intercept array)

%vars
x = 0:numGridsX;                    % X-range
y = 0:numGridsY;                    % Y-range
mb = [m b];                         % Matrix of [slope intercept] values

%functions
lxmb = @(x,mb) mb(1).*x + mb(2);    % Line equation: y = m*x+b
hix = @(y,mb) [(y-mb(2))./mb(1);  y];   % Calculate horizontal intercepts
vix = @(x,mb) [x;  lxmb(x,mb)];    % Calculate vertical intercepts

%calculations for ray going up
L1 = lxmb(x,mb);                    % Calculate Line #1 = y(x,m,b)
hrz = hix(y,mb)';           % [X Y] Matrix of horizontal intercepts
vrt = vix(x,mb)';             % [X Y] Matrix of vertical intercepts
hvix = [hrz; vrt];                 % Concatanated ‘hrz’ and ‘vrt’ arrays

%repeat for the ray going down
m=-numGridsY/(endGrid-startGrid)*2;           % Slope (or slope array)        
b=numGridsY-m*(startGrid);                       % Intercept (or intercept array)
mb = [m b];                         % Matrix of [slope intercept] values
L1 = lxmb(x,mb);                    % Calculate Line #1 = y(x,m,b)
hrz = hix(y,mb)';           % [X Y] Matrix of horizontal intercepts
vrt = vix(x,mb)';             % [X Y] Matrix of vertical intercepts
hvix = [hvix; hrz; vrt];                 % Concatanated ‘hrz’ and ‘vrt’ arrays

%remove out of bound intercepts
exbd = find( (hvix(:,2) < 0) | (hvix(:,2) > max(y)) ); %remove vertical intercepts that are out of bounds
hvix(exbd,:) = [];
exbd = find( (hvix(:,1) < 0) | (hvix(:,1) > max(x)) ); %remove horizontal intercepts that are out of bounds
hvix(exbd,:) = [];
hvix = uniquetol(hvix,1e-5,'ByRows',true);    % Remove repeats (within tolerance) and sort ascending by ‘x’. 

hvdiff=hvix(2:end,:)-hvix(1:end-1,:);
r=sqrt((hvdiff(:,1).*wGrid).^2+(hvdiff(:,2)*hGrid).^2);

a=hvix;
a(:,2)=abs(a(:,2)-max(y));
midPt=(a(2:end,:)-a(1:end-1,:))./2+a(1:end-1,:);
matInd=[midPt(:,2),midPt(:,1)];
matInd=ceil(matInd);


%make matrix mask
%identify which grids the ray passed through
%weight each grid by the path length
myMat=zeros(length(y)-1,length(x)-1);
for k= 1:size(matInd,1)
    myMat(matInd(k,1),matInd(k,2))=myMat(matInd(k,1),matInd(k,2))+r(k); %account for ray passing through same grid 2x
end


%plot intercepts
if display ~=0
    figure(1)
    imagesc([0.5 numGridsX-0.5],[0.5 numGridsY-0.5],flipud(myMat));
%     imagesc(flipud(myMat), 'AlphaData', .1)                         % Draw grids & plot lines
    set(gca,'YDir','normal')
    hold on
    plot(repmat(x,2,length(x)), [0 length(y)-1])    % Vertical gridlines
    hold on
    plot([0 length(x)-1], repmat(y,2,length(y)))    % Horizontal gridlines
    plot([hvix(1,1) hvix(end,1)], [hvix(1,2) hvix(end,2)]) % Plot line
    scatter(hvix(:,1),hvix(:,2))
    hold off
    axis equal
    


end

end
