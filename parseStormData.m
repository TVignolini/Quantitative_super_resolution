function [parsedData] = parseStormData(rawStorm,startFrame)

% Prepares the raw data from Thunderstorm to be processed by other
% functions

parsedData = cell.empty;
frames = 1:max(rawStorm(:,1));
frames = rot90(num2cell(frames),3);

% Allocate space for the matrices to be created next. The maximum possible
% size is too big and fills up all the RAM, so I'm setting an arbitrary
% maximum to the number of localizations per frame which is 10 times the
% average value (rounded to the next integer)
% EDIT 09/11/2020: THIS IS NOT ENOUGH IF AVERAGE LOCALIZATIONS IS LESS THAN
% ONE. Adding a random +5 because otherwise there isd progblim

% EDIT 09/02/2021: THIS IS STILL NOT ENOUGH IF AVERAGE LOCALIZATIONS ARE
% SUPER LOW BUT THERE IS SOME RANDOM SPIKE AT SOME POINT (happening in gel 
% datasets). Turning the +5 in a +20

xPositions = nan(length(frames),ceil((length(rawStorm)/length(frames))*10)+20);
yPositions = nan(length(frames),ceil((length(rawStorm)/length(frames))*10)+20);
sigma = nan(length(frames),ceil((length(rawStorm)/length(frames))*10)+20);
photons = nan(length(frames),ceil((length(rawStorm)/length(frames))*10)+20);
background = nan(length(frames),ceil((length(rawStorm)/length(frames))*10)+20);
tmpM = 1;
tick = 1;
h = waitbar(0,'Parsing raw data');
for n = 1:size(rawStorm,1)
    waitbar(n/size(rawStorm,1))
    
    for m = tmpM:length(frames)
        if isequal(rawStorm(n,1),m)
            if isequal(tmpM,m)
                xPositions(m,tick) = rawStorm(n,2);
                yPositions(m,tick) = rawStorm(n,3);
                sigma(m,tick) = rawStorm(n,4);
                photons(m,tick) = rawStorm(n,5);
                background(m,tick) = rawStorm(n,6);
                tmpM = m;
                tick = tick+1;
            else 
                tick = 1;
                tmpM = m;
                xPositions(m,tick) = rawStorm(n,2);
                yPositions(m,tick) = rawStorm(n,3);
                sigma(m,tick) = rawStorm(n,4);
                photons(m,tick) = rawStorm(n,5);
                background(m,tick) = rawStorm(n,6);
                tick = tick+1;
            end
        end
    end
end

% Remove excess columns
xPositions(:,~any(~isnan(xPositions),1))=[];
yPositions(:,~any(~isnan(yPositions),1))=[];
sigma(:,~any(~isnan(sigma),1))=[];
photons(:,~any(~isnan(photons),1))=[];
background(:,~any(~isnan(background),1))=[];

% Cut the data before the desired starting frame
xPositions(1:startFrame-1,:)=[];
yPositions(1:startFrame-1,:)=[];
sigma(1:startFrame-1,:)=[];
photons(1:startFrame-1,:)=[];
background(1:startFrame-1,:)=[];

parsedData{1} = xPositions;
parsedData{2} = yPositions;
parsedData{3} = sigma;
parsedData{4} = photons;
parsedData{5} = background;
close(h)
end