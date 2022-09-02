%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SHARP seeding algorithm %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% By: Alexander Ernst
% Created: 20Sep20
% Updated: 30Aug22

% OUTPUTS
% isletColor: arbitrarily chosen islet color for COMSOL rendering
% isletDiameters: diameters array, sorted largest to smallest
% isletVolumes: volumes array, sorted largest to smallest
% VMaxO2: maximum oxygen uptake rate
% xPos: x positions, corresponding to isletDiameters and isletVolumes
% yPos: y positions, corresponding to isletDiameters and isletVolumes
% zPos: z positions, corresponding to isletDiameters and isletVolumes

% INPUTS
% cellType: islet source/origin
% denisty: fraction representing volume of device occupied by cells
% roAlg: outer radius of device
% riAlg: inner radius of device
% zAlg: device length

function [isletColor, isletDiameters, isletVolumes, VMaxO2, xPos, yPos, zPos] = sharpSeedingAlgo(cellType, density, roAlg, riAlg, zAlg)

%% 1. Parameters
m_to_um = 1E-6;                                % m to um converter    
m_to_mm = 1E-3;                                % m to mm converter
dMin = 25*m_to_um;                             % min islet diameter (can be lower at the expense of computational burden)
dMax = min(zAlg, roAlg-riAlg)-10*m_to_um;      % max islet diameter (can also be a custom input)
buffer = 15*m_to_um;                           % buffer between islets and walls    
volIEQtgt = density*pi*(roAlg^2-riAlg^2)*zAlg; % target IEQ volume
minTau = 5*m_to_um;                            % minimum viable thickness
fnInd = 0;                                     % master function indexer


%% 2. Run function until it succeeds
while fnInd < 1  % loop until successfully complete
    
    %% 2.1. If some critical thickness tau (500 um) is not met, end simulation                                     
    if roAlg-riAlg<minTau                                               % if actual tau is less than critical tau   
        error('Error. Alginate thickness not sufficient or illogical.') % print this error message
    end                                                                 % end if statement

    
    %% 2.2. Select islet distribution and OCR based on cellType
    switch cellType                                                      % form distribution for cell type entered in input 1 
        case 1                                                           % human islets (cell_type = 1)
            a_HI = 0.36;                                                 % Lognormal shape distribution factor, human islets                                       
            b_HI = 115.4;                                                % Lognormal scale distribution factor, human islets                                    
            PDd = makedist('Lognormal', 'mu', log(b_HI), 'sigma', a_HI); % Human islet distribution
            VMaxO2 = 0.0134;                                             % Human islet OCR (mM/s)
            isletColor = [45/255 194/255 189/255];                       % human islet color
        case 2                                                           % NN islets (cell_type = 2)     
            a_NN = 157.6;                                                % Weibull scale distribution factor, NN islets
            b_NN = 7.78;                                                 % Weibull shape distribution factor, NN islets
            PDd = makedist('Weibull', 'a', a_NN, 'b', b_NN);             % NN islet distribution
            VMaxO2 = 0.0222;                                             % NN islet OCR (mM/s)     
            isletColor = [31/255 79/255 145/255];                        % NN SCB color
        case 3                                                           % MM islets (cell_type = 3)
            a_MM = 254.7;                                                % Weibull scale distribution factor, MM islets
            b_MM = 8.14;                                                 % Weibull shape distribution factor, MM islets
            PDd = makedist('Weibull', 'a', a_MM, 'b', b_MM);             % MM islet distribution
            VMaxO2 = 0.0172;                                             % MM islet OCR (mM/s) 
            isletColor = [8/255 68/255 18/255];                          % MM SCB color
        case 4                                                           % juvenile porcine islets (cell_type = 4)        
            a_PI = 0.36;                                                 % Weibull scale distribution factor, juvenile porcine islets
            b_PI = 152.2;                                                % Weibull shape distribution factor, juvenile porcine islets
            PDd = makedist('Lognormal', 'mu', log(b_PI), 's', a_PI);     % Juvenile porcine islet distribution
            VMaxO2 = 0.0172;                                             % Juvenile porcine islet islet OCR (mM/s) 
            isletColor = [220/255 87/255 31/255];                        % juvenile porcine islet color
        case 5                                                           % rat islets (cell_type = 5)       
            mu_RI = log(112.6);                                          % Lognormal scale distribution factor, rat islets
            sigma_RI = 0.4;                                              % Lognormal shape distribution factor, mouse islets
            PDd = makedist('Lognormal', 'mu', mu_RI, 'sigma', sigma_RI); % Rat islet distribution
            VMaxO2 = 0.034;                                              % Rat islet OCR (mM/s)  
            isletColor = [147/255 40/255 20/255];                        % rat islet color
        case 6                                                           % mouse islets (cellType = 6)      
            a_MI = 168.5;                                                % Weibull scale distribution factor, mouse islets
            b_MI = 3.67;                                                 % Weibull shape distribution factor, mouse islets
            PDd = makedist('Weibull', 'a', a_MI, 'b', b_MI);             % Mouse islet distribution
            VMaxO2 = 0.03;                                               % Mouse islet OCR (mM/s) 
            isletColor = [206/255 202/255 68/255];                       % mouse islet color
        otherwise                                                        % warning message if invalid number is selected      
            error('Error. Invalid cell type selected.');                 % Print error message
    end                                                                  % end switch case cell type distribution selection

    
    %% 2.3. Build a isletDiameters and isletVolumes arrays
    % 2.3.1. Preallocations and loop specific definitions
    volumesPre = zeros(size(linspace(0,1,10000)));            % preallocation 
    diametersPre = zeros(size(linspace(0,1,10000)));          % preallocation
    totVol = 0;                                               % total volume, initialize
    i = 1;                                                    % initialize counter
    % 2.3.2 Create array
    while totVol <= volIEQtgt                                 % loop until sufficient islet volume reached     
        diameterTrial = m_to_um*random(PDd);                  % select islet diameter randomly from PDd
        if (diameterTrial < dMin) || (diameterTrial > dMax)   % if the trial islet is too big or too small...
            continue                                          % ... weed it out         
        end                                                   % end if statement
        diametersPre(i) = diameterTrial;                      % set array of random diameters based on Weibull distribution (each loop picks one number from distribution)
        volumesPre(i) = 4/3*pi*(diametersPre(i)/2)^3;         % calculate volumes array based on diameters (above)                  
        totVol = sum(volumesPre);                             % sum of total volumes
        i = i+1;                                              % increase loop number
    end                                                       % end while loop
    % 2.3.3. Sort arrays
    isletDiameters = sort(nonzeros(diametersPre),'descend')'; % sort diameters in descending order
    isletVolumes = 4/3*pi*(isletDiameters/2).^3;              % convert diamaters to volumes

    
    %% 2.4. Build xPos, yPos, zPos arrays
    % Preallocations and loop-specific definitions
    ONO = [1 -1];                                                                                         % one, negative one (ONO)
    extCoeff = 1.3;                                                                                       % coefficient for increasing the hypothetical seeding box
    j = 1;                                                                                                % indexer
    % % Build x, y, z position arrays, preventing overlap and pokethrough
    while  j <= length(isletDiameters)                                                                    % begin seeding loop
        % Make trial coordinates
        rMaxTrial = roAlg-isletDiameters(j)/2-buffer;                                                     % define maximum r value for this trial islet
        xTrial = ONO(ceil(2*rand))*(rand*(rMaxTrial*extCoeff));                                           % create trial x coordinate, xTrial
        yTrial = ONO(ceil(2*rand))*(rand*(rMaxTrial*extCoeff));                                           % create trial y coordinate, yTrial
        zTrial = buffer+isletDiameters(j)/2+rand*(zAlg-2*(isletDiameters(j)/2+buffer));                   % create trial z coordinate, zTrial
        rTrial = sqrt(xTrial^2+yTrial^2);                                                                 % create representative rTrial
        % Rule out islet seeding outside of the defined dimensions
        if riAlg == 0                                                                                     % if riAlg is 0...                                                                      
            if rTrial >= (roAlg-buffer-isletDiameters(j)/2)                                               % ...and if islet position is greater than ro...
                continue                                                                                  % restart trial seeding process 
            end                                                                                           % end selection if statement      
        else                                                                                              % if riAlg>0...
            if rTrial > (roAlg-buffer-isletDiameters(j)/2) || rTrial < (riAlg+buffer+isletDiameters(j)/2) % ...and if islet position overlaps with ri or ro...
                continue                                                                                  % restart trial seeding process     
            end                                                                                           % end selection if statement
        end                                                                                               % end elimination if statement
        % Check for overlaps
        if j > 2                                                                                          % if seeding the second islet or greater...
            distances = sqrt((xTrial-xPos).^2+(yTrial-yPos).^2+(zTrial-zPos).^2);                         % ...calculate distance between the center of the current and every other previously deposited sphere       
            maxAllowableDistance = (buffer+isletDiameters(j)/2)+isletDiameters(1:(j-1))/2;                % ...calculate array of maximum allowable distances (i.e. the radius of the trial islet and the radius of every other islet plus a buffer)
            differences = distances-maxAllowableDistance;                                                 % ...calculate the differences between the distances between islets and the max allowable distance between islets
            if min(differences) < 0                                                                       % if above array yiels a negative number...
                continue                                                                                  % restart trial seeding process  
            end                                                                                           % end elimination if statement
        end                                                                                               % end seeding if statement
        % Encode islet positions
        xPos(j) = xTrial;                                                                                 % convert successful trial r and phi values to x coordinate
        yPos(j) = yTrial;                                                                                 % convert successful trial r and phi values to y coordinate
        zPos(j) = zTrial;                                                                                 % convert successful trial z value to z coordinate
        j = j+1;                                                                                          % advance loop
    end                                                                                                   % end seeding while loop

    %% 2.5. Double check that no islets overlapped
    % 2.5.1. Preallocations and loop-specific definitions
    mad = zeros(length(isletDiameters),length(isletDiameters));                                 % initialize matrix with length and height of diameters describing the maximum allowable distance
    distances = zeros(length(isletDiameters),length(isletDiameters));                           % initialize matrix with length and height of diameters describing the actual distances
    a = 1;                                                                                      % indexer
    b = 1;                                                                                      % indexer 
    % 2.5.2. Check for overlaps
    for a = 1:length(isletDiameters)                                                            % loop over all diameters in column direction
        for b = 1:length(isletDiameters)                                                        % loop over all diameters in row direction
            if a == b                                                                           % if a = b...
                mad(a,b) = 0;                                                                   % ... the mad is 0 (so put a 0 in the diagonal)
            else                                                                                % else, if not on the diagonal...
                mad(a,b) = (isletDiameters(a)+isletDiameters(b))/2+buffer;                      % ... mad equals the two radii plus a buffer    
            end                                                                                 % end if statement         
            distances(a,b) = sqrt((xPos(a)-xPos(b))^2+(yPos(a)-yPos(b))^2+(zPos(a)-zPos(b))^2); % calculate the actual distance between a and b                               
        end                                                                                     % end row for loop    
        checkRadialDist(a) = roAlg-sqrt(xPos(a)^2+yPos(a)^2)-isletDiameters(a)/2-buffer;        % distance of each islet edge to the alginate edge (including buffer)                                       
    end                                                                                         % end column for loop
    differences = distances - mad;                                                              % calculate difference between actual distance and mad
    if min(min(differences)) < 0 || min(checkRadialDist) < 0                                    % if seeding conditions not met
        fnInd = 0;                                                                              % do not advance loop
        clearvars xPos yPos zPos isletDiameters isletVolumes                                    % clear output variables
    else                                                                                        % if seeding conditions met
        fnInd = 1;                                                                              % advance loop
    end                                                                                         % end if statement
end % end master loop
end % end function