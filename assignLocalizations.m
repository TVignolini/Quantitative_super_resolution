function [localizationsPerCell] = assignLocalizations(numberedCells,coordinates,pxSize,normalizeFlag)

%numberedCells is the output of countBacteria
%coordinates is the output of parseStormData
%pxSize is typically 9.13 (nm) because the masks are enlarged 10x
%normalizeFlag determines if the number of localizations per bacterium is
%going to be divided by the cell area (obtaining a localization density)

cellNumber = max(numberedCells(:));
localizationsPerCell = zeros(cellNumber,1);
%need to create a matrix of localizations by approximating x and y
%positions to pixel number of mask

%extract coordinates from parsed ThunderStorm matrix
xCoordinates = coordinates{1}(:);
xCoordinates = xCoordinates(xCoordinates>0);
yCoordinates = coordinates{2}(:);
yCoordinates = yCoordinates(yCoordinates>0);

for m = 1:length(xCoordinates)
        
    % Convert the coordinates from nm to pixels
    % WARNING: x and y are inverted
    x = yCoordinates(m)/pxSize;
    y = xCoordinates(m)/pxSize;
    approxX = ceil(x);
    approxY = ceil(y);
    if numberedCells(approxX,approxY) ~= 0
        % take the bacterium number in which the localization falls, and
        % add 1 to the number of localization for that bacterium
        localizationsPerCell(numberedCells(approxX,approxY)) = localizationsPerCell(numberedCells(approxX,approxY))+1;
    end
end

if normalizeFlag
    % normalize the number of localizations per bacterium by the surface
    % area of said bacterium
    for m = 1:cellNumber
        [row,col] = find(numberedCells==m,1);
        bacterium = grayconnected(numberedCells,row,col,0);
        sizeBacterium = length(numberedCells(bacterium));
        localizationsPerCell(m) = localizationsPerCell(m)/sizeBacterium;
    end
end
        

end