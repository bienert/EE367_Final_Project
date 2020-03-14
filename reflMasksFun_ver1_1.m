%Nicole Bienert
%Date: 2/24/20

%Purpose: To create reflection coeffecient masks
% glacier surface
% (1,20)  (2,20)  (3,20)  ...  (numGridsX,20)
% (1,19)  (2,19)  (3,19)  ...  (numGridsX,20)
% (1,18)    ...     ...   ...     ...

%version history
%ver1
    %ver1_1: make a matriz of zeros and put a 1 and the bed reflection pt.
    
function myMat = reflMasksFun(txLoc,rxLoc,wGrid,hGrid,numGridsY,numGridsX,display)

startPt=min(txLoc,rxLoc);
endPt=max(txLoc,rxLoc);

endGrid=endPt/wGrid; %number of grids for this antenna separation
startGrid=startPt/wGrid; %Transmitter location normalized to grid width

midPt=((endGrid-startGrid)/2+startGrid);
%make matrix mask
myMat=zeros(numGridsY,numGridsX);
if ceil(midPt)==midPt %on boundary between two grids, weight each grid
    myMat(end,ceil(midPt))=50;
    myMat(end,ceil(midPt)+1)=50;
else
    myMat(end,ceil(midPt))=100;
end

end