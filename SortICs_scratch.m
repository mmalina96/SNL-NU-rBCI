centralICs = zeros(20,128);
for i = 1:20
    for j = 1:128
        if centralLocs(i,j) == 0
            centralICs(i,j) = 1;
        elseif centralLocs(i,j) < 50
            val = (-centralLocs(i,j)^2)+1;
            centralICs(i,j) = val / 2500;
        else
            centralICs(i,j) = 0;
        end
    end
end

    

d = -50:.0001:50;
y = ((-d.^2)+1) ./ 2500;
plot(d,y)

neighborVec = zeros(20,128);
for i = 1:128
    for j = 1:20
         if centralLocs(j,i) < 50
             neighborVec(j,i) = CNindex(j,i);
         else
             neighborVec(j,i) = 0;
         end
    end
end

centralICs = zeros(20,128);
for i = 1:20
    for j = 1:128
        if centralLocs(i,j) == 0
            centralICs(i,j) = 1;
        elseif centralLocs(i,j) < 50
            val = (centralLocs(i,j)^2)+1;
            centralICs(i,neighborVec) = val / 2500;
        else
            centralICs(i,j) = 0;
        end
    end
end

euDistNorm = zeros(128,128);
for i = 1:128
    euDistNorm(i,:) = (distVec(i,i) - min(distVec(i,:)) ./ max(distVec(i,:)) - min(distVec(i,:)));
end

euDistNorm = normalize(distVec, 2, 'range', [0,1]);

%% fixing the iterative part




centralICs = zeros(20,128);
d = [1:20:901];
for k = 1:size(d,2)-1
    for i = 1:20
        for j = 1:128
            if centralLocs(i,j) == 0
                centralICs(i,j) = 1;
            elseif  centralLocs(i,j) < 0.3 % what if you came up with a vector 
                    % that was like (0.1:0.1:0.5) to rand determine the how far
                    % out the effected electrodes go. Might make generalization
                    % better
                    centralICs(i,j) = (-(centralLocs(i,j)^(0.5))+1); % dist 1
                   % centralICs(i,j) = (-(centralLocs(i,j)^(0.25))+1); % dist 2
                   % centralICs(i,j) = (-(centralLocs(i,j)^(0.1)+1); % dist 3
                   % centralICs(i,j) = (-(centralLocs(i,j) + 1); % dist 4
                   % centralICs(i,j) = exp(-(centralLocs(i,j)).^2) % normal 
                   % centralICs(i,j) = exp(-(centralLocs(i,j)).^4) % normal 
            else
                centralICs(i,j) = 0;
            end
        end
    end
    keepMat(d(k):(d(k+1)-1),:) = centralICs;    
end
































