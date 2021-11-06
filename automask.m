function [mask, locDensity, locError] = automask(coordinates,maskSizeX,maskSizeY,pxSize,radius,lothreshold,hithreshold,analyzeFlag,offset,refmask)

%Creates masks based on the localization maps. Creates a round mask, 1 
%micron (can be defined) in diameter centered around every localization.
%Rationale: if a bacterium is labelled, on focus, and alive at the time of 
%fixation, it should contain at least one visible PAmCherry. Excludes empty
%areas and false negatives are minimized (fluorescent areas overlap well
%with on-focus areas).

%Mk2 edit: every circle now drawn onto the map adds to the value of each
%pixel it intersects. A treshold can be implemented to remove areas below a
%specific density. Excludes too much sparsely-populated areas which could
%be over-represented otherwise.

%Mk3 edit: this function can now also give the final localization density
%in the resulting mask area, without the need to process the data through
%all the other functions. If analyzeFlag = 1 the output includes locDensity
%in the form of localizations per square micron

%Mk4 edit: upper threshold inserted, to remove torches

%Mk5 edit: offset inserted to remove background localizations (if blank
%calibration was performed). If no offset is known, input 0

%Mk6 edit: refmask inserted to avoid redrawing mask every time if one is
%already present in the workspace. If no mask is present, input 0

%Mk7 edit (02/02/2021): added a Poissonian error estimation based on number
%of localizations

%coordinates is the output of parseStormData
%maskSizeX and Y are typically either 512 or 5120 (original image size is
%512x512 px)
%pxSize should be 9.13 (nm) if the masks are 5120x5120 or 91.3 if the
%masks are 512x512

radius = radius/pxSize;

%extract coordinates from parsed ThunderStorm matrix
xCoordinates = coordinates{1}(:);
xCoordinates = xCoordinates(xCoordinates>0);
yCoordinates = coordinates{2}(:);
yCoordinates = yCoordinates(yCoordinates>0);


if ~any(refmask(:))
mask = zeros(maskSizeY,maskSizeX,'uint16');
[maskColumns, maskRows] = meshgrid(1:maskSizeX, 1:maskSizeY);

%draw circular mask of defined radius around each localized spot and add
%one on top of another sequentially

h = waitbar(0,'Drawing circles');
for m = 1:length(xCoordinates)
    waitbar(m/length(xCoordinates))
    x = xCoordinates(m)/pxSize;
    y = yCoordinates(m)/pxSize;
    circle = (maskRows - y).^2 + (maskColumns - x).^2 <= radius.^2;
    circle = im2uint16(circle)/65535;
    mask = mask+circle;
end

%prune the mask where the localization density is too low or too high
masklo = mask>lothreshold;
maskhi = mask<hithreshold;
mask = masklo.*maskhi;

%convert uint16 mask into uint8 for compatibility with assignlocalizations
mask = im2uint8(logical(mask))/255;

close(h)

else
    mask = refmask;
end

if analyzeFlag
    
    %convert localizations per pixel to localizations per square micron
    for m = 1:length(xCoordinates)
        
    % Convert the coordinates from nm to pixels
    % WARNING: x and y are inverted
    x = xCoordinates(m)/pxSize;
    y = yCoordinates(m)/pxSize;
    approxX = ceil(x);
    approxY = ceil(y);
    if mask(approxY,approxX) == 0
        % take the bacterium number in which the localization falls, and
        % add 1 to the number of localization for that bacterium
        xCoordinates(m) = 0;
    end
    end
    
    xCoordinates = xCoordinates(xCoordinates>0);
    
    maskArea = ((pxSize^2)*nnz(mask))/1000000;
    locDensity = (length(xCoordinates)-offset)/maskArea;
    % associate Poissonian error to number of localizations
    locError = sqrt(length(xCoordinates)-offset)/maskArea;
else
    locDensity = "None";
    locError = "None";
end

end
    