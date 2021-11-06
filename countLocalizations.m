function [locsPerFrame] = countLocalizations(mask,coordinates,pxSize)
%Creates a vector of localizations per frame, considering only
%localizations happening within the given mask

%coordinates is the output of parseStormData
%pxSize is typically 9.13 (nm) because the masks are enlarged 10x

locsPerFrame = zeros(size(coordinates{1},1),1);

for n = 1:length(locsPerFrame)
    for m = 1:length(coordinates{1}(n,~isnan(coordinates{1}(n,:))))
        
        % Convert the coordinates from nm to pixels
        x = coordinates{2}(n,m)/pxSize;
        y = coordinates{1}(n,m)/pxSize;
        approxX = ceil(x);
        approxY = ceil(y);
        if mask(approxX,approxY) ~= 0
            % take the bacterium number in which the localization falls, and
            % add 1 to the number of localization for that bacterium
            locsPerFrame(n) = locsPerFrame(n)+1;
        end
    end
end        

end