function [numberedMask] = countBacteria(mask,maxValue,sizeThreshold)

%sizeThreshold: minimum size for particle to be considered a cell

%Eats up a segmented picture and assigns a different value to each
%separated particle. Only accepts uint8 matrices

b = 1; %bacteria counter

h = waitbar(0,'looking for bugs');

for n = 1:size(mask,2) %columns
    
    waitbar(n/size(mask,2))
    
    for m = 1:size(mask,1) %rows
        if mask(m,n) == maxValue %check if there is something here
            bacterium = grayconnected(mask,m,n,0);
            sizeBacterium = length(mask(bacterium));
            bacterium = uint8(bacterium);
            mask = mask-bacterium; %turn all the 255s in the newly found bacterium into 254s
            for j = 1:size(mask,2)
                for k = 1:size(mask,1)
                    if mask(k,j) == (maxValue-1) %turn all the 254s in b
                        if sizeBacterium > sizeThreshold %check if it's a real bacterium or a small artifact
                            mask(k,j) = b;
                        else
                            mask(k,j) = 0;
                        end
                    end
                end
            end
            if sizeBacterium > sizeThreshold 
                b = b+1;
            end
        end
    end
end

close(h)

numberedMask = mask;

end