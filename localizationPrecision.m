function [Uncertainties,signalComponent,backgroundComponent] = localizationPrecision(localizations,range,threshold,patience,pxSize,neighbourhoodWatch)

% Last update: 27/09/2018

% Newest version: calculate Thompson components on every spot, not only the
% repeated ones.

% "localizations" is the output of parseStormData

% "threshold" is the max allowed distance (in nm) between spots in
% consecutive frames under which they can be considered to be the same
% molecule

%% Identification of persisting fluorophores

xPositions = localizations{1};
yPositions = localizations{2};
sigma = localizations{3};
photons = localizations{4};
background = localizations{5};
frames = 1:size(xPositions,1);
% Prepare empty matrices to be used next
xUncertainties = [];
yUncertainties = [];
Uncertainties = [];
forbiddenCoordinates = [0 0];
forbiddenCoordinatesTmp = [];
signalBucket = [];
backgroundBucket = [];
signalComponent = [];
backgroundComponent = [];
% For every localized spot, look for the nearest localization in the
% following frame. If the distance is lower than the threshold, start from
% the newest localization and move on to the next frame

h = waitbar(0,'Looking for persisting fluorophores');
% Start from the first frame and go on until the last possible frame that
% allows a search window equal to the specified range
for n = 1:(length(frames)-(range-1))
    waitbar(n/(length(frames)-(range-1)))
    % Within a frame, analyze all of the localized spots from first to last
    for m = 1:nnz(~isnan(xPositions(n,:)))
        patienceTmp = patience;
        gap = 0;
        check1 = [];
        check2 = [];
        MinIndex = [];   
        % Check if the selected spot has already been used for an accuracy
        % estimation. If it has, move on to the next localization within 
        % the frame
        check1 = unique(find(forbiddenCoordinates(:,1)==n));
        check2 = unique(find(forbiddenCoordinates(:,2)==m));
        if isequal(length(cat(1,check1,check2)),length(unique(cat(1,check1,check2))))
            check1 = [];
            check2 = [];
            positions = [];
            MinIndex(1,1) = m;
            positions(1,1) = xPositions(n,m);
            positions(1,2) = yPositions(n,m);
            for a = 0:(range-2)
                % Check for gap frames
                if (a > 0 && isequal(positions(a,1),positions(a+1,1)))
                    xDistances = abs(xPositions(n+a-(gap),MinIndex(a+1,1))-xPositions(n+a+1,:));
                    yDistances = abs(yPositions(n+a-(gap),MinIndex(a+1,1))-yPositions(n+a+1,:));
                    Distances = [];
                    for k = 1:length(xDistances)
                        Distances = [Distances;sqrt(xDistances(k)^2+yDistances(k)^2)];
                    end
                    [MinValue,MinIndex(a+2,1)] = min(Distances);
                    MinValue2 = min(setdiff(Distances(:),MinValue));
                    if isempty(MinValue2)
                        MinValue2 = NaN;
                    end
                else
                    % Generates a vector of distances between the selected
                    % spot and all of the spots in the next frame. 
                    xDistances = abs(xPositions(n+a,MinIndex(a+1,1))-xPositions(n+a+1,:));
                    yDistances = abs(yPositions(n+a,MinIndex(a+1,1))-yPositions(n+a+1,:));
                    Distances = [];
                    for k = 1:length(xDistances)
                        Distances = [Distances;sqrt(xDistances(k)^2+yDistances(k)^2)];
                    end
                    [MinValue,MinIndex(a+2,1)] = min(Distances);
                    MinValue2 = min(setdiff(Distances(:),MinValue));
                    if isempty(MinValue2)
                        MinValue2 = NaN;
                    end
                end
                % Same check as last one
                check1 = unique(find(forbiddenCoordinates(:,1)==n+a+1));
                check2 = unique(find(forbiddenCoordinates(:,2)==MinIndex(a+2,1)));
                if ~isequal(length(cat(1,check1,check2)),length(unique(cat(1,check1,check2))))
                    check1 = [];
                    check2 = [];
                    break
                end
                if MinValue > threshold || isnan(MinValue)
                    patienceTmp = patienceTmp-1;
                    gap = gap+1;
                    positions(a+2,1) = positions(a+1,1);
                    positions(a+2,2) = positions(a+1,2);
                    MinIndex(a+2,1) = MinIndex(a+1,1);                                       
                else
                    if neighbourhoodWatch && ((MinValue2 <= threshold))
                        break
                    end
                    positions(a+2,1) = xPositions(n+a+1,MinIndex(a+2,1));
                    positions(a+2,2) = yPositions(n+a+1,MinIndex(a+2,1));
                    tmpRow = n+a+1;
                    tmpColumn = MinIndex(a+2,1);
                    forbiddenCoordinatesTmp = [forbiddenCoordinatesTmp;tmpRow tmpColumn];
                    gap = 0;
                end
                if patienceTmp < 0
                    break
                end                
                if isequal (a,range-2)
                    xUncertainties = [xUncertainties;std(positions(:,1))];
                    yUncertainties = [yUncertainties;std(positions(:,2))];
                    Uncertainties = mean(cat(2,xUncertainties,yUncertainties),2);
                    forbiddenCoordinates = [forbiddenCoordinates;forbiddenCoordinatesTmp];
                    if isnan(std(unique(positions(:,1)))) || isnan(std(unique(positions(:,2))))
                        positions(:,1);
                        positions(:,2);
                        xPositions(n,m)
                        yPositions(n,m)
                    end
                        
                end
            end
        end
        MinIndex = [];        
        forbiddenCoordinatesTmp = [];
    end
end
close(h)
% Calculate Thompson accuracy on localized spots
for c = 1:length(frames)
    for b = 1:nnz(xPositions(c,:))
        signalTmp = (2*(sigma(c,b)^2)+((pxSize^2)/12))/photons(c,b);
        backgroundTmp = (8*pi*(sigma(c,b)^4)*(background(c,b)^2))/((pxSize^2)*(photons(c,b)^2));
        signalBucket = [signalBucket;signalTmp];
        backgroundBucket = [backgroundBucket;backgroundTmp];
    end
end
signalComponent = [mean(signalBucket),std(signalBucket)];
backgroundComponent = [mean(backgroundBucket),std(backgroundBucket)];
end