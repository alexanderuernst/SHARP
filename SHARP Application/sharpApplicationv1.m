%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SHARP %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Author: Alexander Ulrich Ernst, PhD (email: aue2@cornell.edu)
% Date created: 20Jul20
% Last updated: 28Jul22

% This program calculates oxygen distributions and related outcomes in user
% -programmable bioartificial pancreas devices comprising hydrogel and 
% islets.

% INSTRUCTIONS
% 1. Follow input prompts to enter device contents: cell type, cell volume, 
% device type, device parameters, and number of monte carlo iterations you
% wish the program to run. The program will calculate three device fates:
% mean pO2, net loss of insulin secretory capacity (Psi), and net
% volumetric fraction of necrotic tissue of the islets. It will calculate
% these outcomes on a population basis (e.g., the average pO2 of all the
% encapsulated islets) and on an individual islet basis. The former is
% printed as a table and the latter is saved as matlab files to a
% designated folder. In addition, the program will display the device
% geometry for each iteration of the simulation and will save the first
% iteration as a COMSOL file for the user's convenience. 

% 2. If the user wishes to change the values of individual parameters  
% (e.g., the value of the external oxygen concentrioton), this can be done 
% manually in section 5.6. wherein the parameter values are hard-coded. The 
% default path is set to the author's Dropbox, but can (and should) be 
% changed manually by the user.


%% 1. Clear previous structures 
clc
clear all
close all


%% 2. Randomize seed
rng shuffle
rand;


%% 3. Assign application folder (used for future saving functions)
sharpFolder = 'C:\Users\Alexander Ernst\Dropbox\SHARP Application';


%% 4. Parameters
% 4.1. Conversion factors
m2um = 10^-6;               % conversion factor from m to um
m2mm = 10^-3;               % conversion factor from m to mm
m2cm = 10^-2;               % conversion factor from m to cm
% 4.2. General islet information
buffer = 10*m2um;           % buffer between islets
d1IEQ = 150*m2um;           % diameter of standard 150 um islet
v1IEQ = 4/3*pi*(d1IEQ/2)^3; % volume of one standard 150 um diameter islet
% 4.3. Monte Carlo loop preinitializers
ind = 0;                    % initialize master indexer
didFail = 0;                % initialize failure indicator 


%% 5. Inputs
inputFailure = 0;                                                                                                 % initialize input failure indicator
while inputFailure < 1                                                                                            % begin input while loop
    % 5.1. BAP geometry structure and cell type
    bapGeometry = input('Enter BAP geometry type (1 = planar slab, 2 = cylinder, 3 = hollow cylinder):  ');       % geometry of BAP
    cellType = input('Enter islet source (1 = human, 2 = NN SCBs, 3 = MM SCBs, 4 = pig, 5 = rat, 6 = mouse):  '); % cell source
    if rem(bapGeometry,1) == 0 && rem(cellType,1) == 0 && bapGeometry <= 3 && cellType <= 6                       % if inputs are valid...
        inputFailure = 1;                                                                                         % ...advance the input failure indicator
    else                                                                                                          % otherwise...
        inputError1Readout = '\n\nError: invalid input(s), please reenter\n\n\n';                                 % ...readout statement  (to be printed)
        fprintf(inputError1Readout)                                                                               % ...print readout
        continue                                                                                                  % ...return to top of loop
    end                                                                                                           % end if statement
    % 5.2. Islet encapsulated islet payload
    targetIEQ = input('Enter the total islet volume, in IEQ, encapsulated in the BAP device:  ');                 % total islet payload
    % 5.2. Geometry parameters
    switch bapGeometry                                                                                            % begin switch/case to assign parameters for BAP device geometry
        case 1                                                                                                    % case 1 = slab
            zAlg = m2um*input('Enter the thickness of the hydrogel slab (in micrometers):  ');                    % slab thickness (if using primary islets, must be >= 500 um)
            riAlg = 0;                                                                                            % inner radius of slab (= 0)
            roAlg = m2mm/2*input('Enter the diameter of the hydrogel slab (in millimeters):  ');                  % radius of slab
            volAlg = pi*roAlg^2*zAlg;                                                                             % device total volume
            targetDensity = v1IEQ*targetIEQ/volAlg;                                                               % volumetric islet density calculation
        case 2                                                                                                    % case 2 = cylinder (if using primary islets, must be >= 500 um)              
            roAlg = m2um/2*input('Enter the diameter of the hydrogel (in micrometers):  ');                       % cylinder radius (if using primary islets, must be >= 500 um)
            riAlg = 0;                                                                                            % inner radius of cylinder (= 0)
            zAlg = m2mm*input('Enter the length of the hydrogel cylinder (in millimeters):  ');                   % length of cylinder
            volAlg = pi*roAlg^2*zAlg;                                                                             % device total volume
            targetDensity = v1IEQ*targetIEQ/volAlg;                                                               % volumetric islet density calculation
        case 3                                                                                                    % case 3 = hollow cylinder
            riAlg = m2um/2*input('Enter the inner diameter of the hollow cylinder (in micrometers):  ');          % hollow cylinder inner radius
            tauAlg = m2um*input('Enter the thickness of the hydrogel concentric core (in micrometers):  ');       % hollow cylinder hydrogel thickness (if using primary islets, must be >= 500 um)
            roAlg = riAlg+tauAlg;                                                                                 % total hollow cylinder radius
            zAlg = m2mm*input('Enter the length of the hollow cylinder (in millemeters):  ');                     % length of hollow cylinder
            volAlg = pi*(roAlg^2-riAlg^2)*zAlg;                                                                   % device total volume
            targetDensity = v1IEQ*targetIEQ/volAlg;                                                               % volumetric islet density calculation
    end                                                                                                           % end geometry parameter assignment
    densityReadout = ['\nThe cell density is ' num2str(round(100*targetDensity,4)) '%% v/v \n\n'];                % write volumetric islet density readout for confirmation  
    fprintf(densityReadout);                                                                                      % print density readout
    % 5.4. Enter number Monte Carlo iterations (>=3 recommended, or a minimum of 1500 total IEQ simualted)
    nRuns = input('Enter number of iterations you wish the Monte Carlo simulation to perform:  ')';               % # loops in MC simulation
end                                                                                                               % end input while loop


%% 5. Run MC simulation
while ind < nRuns

    %% 5.1. Advance (or don't advance) indexer
    % 5.1.1. Check failure status
    switch didFail                                                      % begin switch/case to handle iteration numbers given failures
        case 0                                                          % if no error occurred...
            ind = ind+1;                                                % ... advance indexer
        case 1                                                          % if  error occured
            ind = ind+0;                                                % do not advance indexer
    end                                                                 % end switch case
    % 5.1.2. Print iteration number
    iteration = ['\n<strong>Iteration #: ' num2str(ind) '</strong>\n']; % create iteration string
    fprintf(iteration)                                                  % print iteration string
    didFail = 0;                                                        % reset failure indicator
    

    %% 5.2. Load COMSOL
    import com.comsol.model.*              % COMSOL class
    import com.comsol.model.util.*         % COMSOL class 
    model = ModelUtil.create('model');     % Create matlab variable model    
    geom1 = model.geom.create('geom1', 3); % Create 3D geometry
    ModelUtil.showProgress(true);          % Display progress bar

    
    %% 5.3. Implement global geometry parameters
    model.param.set('zAlg', [num2str(zAlg) '[m]'], 'length of alginate radius');            % create comsol parameter zAlg
    model.param.set('roAlg', [num2str(roAlg) '[m]'], 'outer radius of the alginate layer'); % create comsol parameter roAlg 
    model.param.set('riAlg', [num2str(riAlg) '[m]'], 'outer radius of the alginate layer'); % create comsol parameter roAlg 
    model.param.set('buffer', [num2str(buffer) '[m]'], 'buffer');                           % create comsol parameter buffer
   
    
    %% 5.4. Build geometry as COMSOL entity
    % 5.4.1. Build hydrogel domain
    alginateDomain = model.geom('geom1').create('cyl1', 'Cylinder');                                                                      % create the aglinate cylinder
    alginateDomain.set('r', 'roAlg');                                                                                                     % set radius
    alginateDomain.set('h', 'zAlg');                                                                                                      % set height
    alginateDomain.set('selresult', true);                                                                                                % allow to appear in selections
    alginateDomain.set('selresultshow', 'all');                                                                                           % add all aspects to selections
    alginateDomain.set('color', 'custom');                                                                                                % Make custom color selection
    alginateDomain.set('customcolor', [234/255 242/255 242/255]);                                                                         % change the displayed color of the cylinder
    % 5.4.2. Build islets
    [isletColor, isletDiameters, isletVolumes, VMaxO2, xPos, yPos, zPos] = sharpSeedingAlgo(cellType, targetDensity, roAlg, riAlg, zAlg); % create array of islets and their positions and OCRs                                                                                                                   % end islet color switch/case
    model.geom('geom1').selection.create('csel1', 'CumulativeSelection');                                                                 % cumulative selection of all islets
    k = 1;                                                                                                                                % indexer
    for k = 1:length(isletDiameters)                                                                                                      % loop over all islets    
        islettag = model.geom('geom1').feature().uniquetag('sph');                                                                        % create a unique name for a sphere i.e. sph(k)      
        islet{k} = model.geom('geom1').feature().create(islettag, 'Sphere');                                                              % create a sphere
        islet{k}.set('r', isletDiameters(k)/2);                                                                                           % set radius of sphere to indexed diameter value
        islet{k}.set('pos', [xPos(k) yPos(k) zPos(k)]);                                                                                   % set position of sphere to indexed coordinate value
        islet{k}.set('contributeto','csel1');                                                                                             % contribute the new islet to the cumulative selection that is all islets
        islet{k}.set('selresult', true);                                                                                                  % "  
        islet{k}.set('selresultshow', 'all');                                                                                             % add all islet features to selections  
        islet{k}.set('color','custom');                                                                                                   % change the displayed color of the islets
        islet{k}.set('customcolor', isletColor);                                                                                          % set islet color for rendering in COMSOL
    end                                                                                                                                   % end loop
    % 5.4.3. Print geometry
    model.geom('geom1').run;                                                                                                              % run the geometry
    model.geom('geom1').run('fin');                                                                                                       % form union
    figure(1)                                                                                                                             % open new figure window
    mphgeom(model,'geom1', 'facealpha', .5)                                                                                               % Plot geometry, make it sort of transparent
    
    
    %% 5.5. Create face selections
    % 5.5.1. Cylinder side faces
    model.geom('geom1').selection().create('csel2', 'CumulativeSelection');               % create cumulative selection 2 for cylinder sides
    ballsel1 = model.component('mod1').geom('geom1').create('ballsel1', 'BallSelection'); % create ball selection subnode
    ballsel1.set('entitydim', 2);                                                         % select faces  
    ballsel1.set('r', 'buffer/2');                                                        % select size of ball selector
    ballsel1.set('posz', 'zAlg/2');                                                       % set z position of ball selector
    ballsel1.set('posy', '0');                                                            % set x position of ball selector
    ballsel1.set('posx', 'roAlg');                                                        % set x position of ball selector
    ballsel1.set('condition', 'intersects');                                              % create selection for all entities within ball
    ballsel1.set('contributeto', 'csel2');                                                % contribute to cumulative selection 2
    ballsel2 = model.component('mod1').geom('geom1').create('ballsel2', 'BallSelection'); % create ball selection subnode
    ballsel2.set('entitydim', 2);                                                         % select faces  
    ballsel2.set('r', 'buffer/2');                                                        % select size of ball selector
    ballsel2.set('posz', 'zAlg/2');                                                       % set z position of ball selector
    ballsel2.set('posy', '0');                                                            % set x position of ball selector
    ballsel2.set('posx', '-roAlg');                                                       % set x position of ball selector
    ballsel2.set('condition', 'intersects');                                              % create selection for all entities within ball
    ballsel2.set('contributeto', 'csel2');                                                % contribute to cumulative selection 2
    % 5.5.2. Cylinder end faces
    model.geom('geom1').selection().create('csel3', 'CumulativeSelection');               % create cumulative selection 2 for cylinder sides
    ballsel3 = model.component('mod1').geom('geom1').create('ballsel3', 'BallSelection'); % create ball selection subnode
    ballsel3.set('entitydim', 2);                                                         % select faces  
    ballsel3.set('r', 'buffer/2');                                                        % select size of ball selector
    ballsel3.set('posz', '0');                                                            % set z position of ball selector
    ballsel3.set('posy', 'riAlg+(roAlg-riAlg)/2');                                        % set x position of ball selector
    ballsel3.set('posx', '0');                                                            % set x position of ball selector
    ballsel3.set('condition', 'intersects');                                              % create selection for all entities within ball  
    ballsel3.set('contributeto', 'csel3');                                                % contribute to cumulative selection 3
    ballsel4 = model.component('mod1').geom('geom1').create('ballsel4', 'BallSelection'); % create ball selection subnode
    ballsel4.set('entitydim', 2);                                                         % select faces  
    ballsel4.set('r', 'buffer/2');                                                        % select size of ball selector
    ballsel4.set('posz', 'zAlg');                                                         % set z position of ball selector
    ballsel4.set('posy', 'riAlg+(roAlg-riAlg)/2');                                        % set x position of ball selector
    ballsel4.set('posx', '0');                                                            % set x position of ball selector
    ballsel4.set('condition', 'intersects');                                              % create selection for all entities within ball 
    ballsel4.set('contributeto', 'csel3');                                                % contribute to cumulative selection 3

    
    %% 5.6. Implement variables
    % 5.6.1. Variables that live in all domains
    allVars = model.component('mod1').variable.create('var1');                                                    % create variable node
    allVars.selection.named('geom1_cyl1_dom');                                                                    % assign whole geometry to node            
    allVars.set('D_o2_alg', '2.7E-9 [m^2/s]', 'diffusivity of oxygen in alginate');                               % alginate diffusivity
    allVars.set('alpha_o2_alg', '9.3E-6 [mol/m^3/Pa]', 'Bunsen solubility of oxygen in alginate');                % alginate o2 solubility
    allVars.set('Pext', '40 [mmHg]', 'external pO2');                                                             % external pO2
    allVars.set('pO2', 'c/alpha_o2_alg', 'conversion to partial pressure');                                       % pO2 definition
    % 5.6.2. Variables that live in the islets
    isletVars = model.component('mod1').variable.create('var2');                                                  % create variable node
    isletVars.selection.named('geom1_csel1_dom');                                                                 % assign islets to node
    isletVars.set('D_o2_islet', '2.0E-9 [m^2/s]', 'diffusivity of oxygen in islets');                             % islet diffusivity
    isletVars.set('Vmax', [num2str(VMaxO2) '[mM/s]'], 'maximum oxygen uptake in islets');                         % islet max OCR
    isletVars.set('Km', '1E-3 [mM]', 'half-maximal coefficient');                                                 % half-maximal constant
    isletVars.set('cCrit', '1E-4 [mM]', 'critical oxygen concentration required for survival');                   % critical oxygen level
    isletVars.set('Vo2', '-Vmax*c/(c+Km)*flc1hs((c-cCrit)[1/mM],5E-5)', 'Real OCR');                              % islet OCR (Michaelis-Menten kinetics)
    isletVars.set('KmIns', '2 [mmHg]', 'half-max. coeff., oxygen-modulated insulin release');                     % Km, O2-Ins
    isletVars.set('sIns', '(pO2 [1/mmHg])^3/((pO2 [1/mmHg])^3+(KmIns [1/mmHg])^3)*flc1hs((c-cCrit)[1/mM],5E-5)'); % insulin secretion second phase potential

    
    %% 5.7. Create nonuniform mesh
    mesh1 = model.mesh.create('mesh1','geom1');              % create a mesh node
    ftet1 = mesh1.feature().create('ftet1', 'FreeTet');      % create a free tetrahydral meshing process
    ftet1.selection.geom('geom1', 3);                        % select a geometry of which to mesh in ftet1
    ftet1.selection.all;                                     % select all islets to mesh
    size1 = ftet1.create('size1', 'Size');                   % create a size subnode of the mesh to tune mesh building parameters individually
    size1.selection.all;                                     % select all the islets to mesh
    size1.set('custom', 'on');                               % build a custom mesh
    size1.set('hmax', 250E-5);                               % set maximum element size (m)
    size1.set('hmaxactive', true);                           % active above value    
    size1.set('hmin', 4.0E-7);                               % set minimum element size (m)
    size1.set('hminactive', true);                           % active above value
    size1.set('hcurve', 0.26);                               % set curvature factor
    size1.set('hcurveactive', true);                         % activate above value
    size1.set('hnarrow', 2.5);                               % set resolutoin of narrow regions value
    size1.set('hnarrowactive', true);                        % activate above value
    size1.set('hgrad', 1.32);                                % set maximum element growth rate
    size1.set('hgradactive', true);                          % activate above value
    try                                                      % try/catch statement to catch meshing errors   
        model.mesh('mesh1').run;                             % run the mesh   
    catch                                                    % If there's an error...    
        didFail = 1;                                         % 1 indicates that the meshing did fail 
        meshFailureMessage = 'Failed to mesh.\n';            % seeding failure readout
        fprintf(meshFailureMessage)                          % print readout
        clearvars xPos yPos zPos isletDiameters isletVolumes % ...clear all the islet variables including the master diameters list
        continue                                             % return to top of loop
    end                                                      % end meshing trycatch

    
    %% 5.8. Print specs
    % 5.8.1. Salient calculations & recordings
    ieqIteration = sum(isletVolumes)/v1IEQ;                                     % IEQ of islets
    allISI(ind) = ieqIteration/length(isletDiameters);                          % ISI of iteration
    actualDensity = sum(isletVolumes)/volAlg;                                   % actual density     
    densityIteration(ind) = targetDensity*100;                                  % track the input density each iteration (%)
    % 5.8.2. Write and print specs
    inputSpecs = ['Mesh succesful' '\n'...                                      % indicate mesh succesfully completed
                  'ISI: ' num2str(allISI(ind)) '\n'...                          % ISI of this iteration
                  'Target IEQ: ' num2str(targetIEQ) '\n'...                     % target IEQ
                  'Actual IEQ: ' num2str(ieqIteration) '\n'...                  % actual IEQ
                  'Target density: ' num2str(targetDensity*100) '\n'...         % target density
                  'Actual density: ' num2str(actualDensity*100) '\n'...         % actual density
                  'Max diameter: ' num2str(max(isletDiameters)/m2um) '\n'...    % max diameter
                  'Min diameter: ' num2str(min(isletDiameters)/m2um) '\n'...    % min diameter
                  'Number of islets: ' num2str(length(isletDiameters)) '\n\n']; % number of islets
    fprintf(inputSpecs)                                                         % print spec string

    
    %% 5.9. Implement oxygen physics
    % 5.9.1. TDS master subnode
    tds = model.physics.create('tds', 'DilutedSpecies', 'geom1');                                             % create transport of diluted species (TDS) physics node
    tds.prop('TransportMechanism').set('Convection', false);                                                  % ... omit convection
    % 5.9.2. CDM subnode - alginate
    tds.feature('cdm1').set('D_c', {'D_o2_alg'; '0'; '0'; '0'; 'D_o2_alg'; '0'; '0'; '0'; 'D_o2_alg'});       % set the diffusivity of alginate    
    % 5.9.3. CDM subnode - islets
    tds.create('cdm2', 'ConvectionDiffusionMigration', 3);                                                    % create a CDM subnode representing the islets
    tds.feature('cdm2').selection.named('geom1_csel1_dom');                                                   % assign the islets to the CDM2 subnode 
    tds.feature('cdm2').set('D_c', {'D_o2_islet'; '0'; '0'; '0'; 'D_o2_islet'; '0'; '0'; '0'; 'D_o2_islet'}); % set the diffusivity of islet tissue 
    % 5.9.4. Initial conditions subnode
    tds.feature('init1').set('initc', 'Pext*alpha_o2_alg');                                                   % set the initial concentration in alginate domains to a value > 0
    initc2 = tds.feature().create('initc2', 'init',3);                                                        % create new initial condition node for islets 
    initc2.selection.named('geom1_csel1_dom');                                                                % select the islets (csel1)  
    initc2.set('initc', 'cCrit');                                                                             % set initial value to c_critical 
    % 5.9.5. RXN subnode - islet oxygen consumption
    tds.create('reac1', 'Reactions', 3);                                                                      % create a reaction term subnode
    tds.feature('reac1').selection.named('geom1_csel1_dom');                                                  % assign the islets to the reaction subnode
    tds.feature('reac1').set('R_c', 'Vo2');                                                                   % set the reaction term to the MM kinetics multiplied by a heaviside function    
    % 5.9.6. Concentration - host interface
    peritonealInterface = tds.create('conc1', 'Concentration', 2);                                            % create a boundary condition subnode for the outer boundary                                                              
    switch bapGeometry                                                                                        % switch/case to assign host interface depending on bap geometry
        case 1                                                                                                % case 1 - planar slab
            peritonealInterface.selection.named('geom1_csel3_bnd');                                           % assign the hydrogel slab end faces to the boundary condition subnode
        case 2                                                                                                % case 2 - cylinder
            peritonealInterface.selection.named('geom1_csel2_bnd');                                           % assign hydrogel cylinder curved faces to the boundary condition subnode
        case 3                                                                                                % case 3 - hollow cylinder
            peritonealInterface.selection.named('geom1_csel2_bnd');                                           % assign hydrogel cylinder cruved faces to the boundary condition subnode
    end                                                                                                       % end switch/case
    peritonealInterface.set('species', true);                                                                 % apply the boundary condition to the surface
    peritonealInterface.set('c0', 'Pext*alpha_o2_alg');                                                       % set the boundary condition to the product of the solubility and partial pressure
    % 5.9.7. No flux subnode
    noFluxEnds = tds.create('nflx', 'NoFlux', 2);                                                             % create a surface reaction term subnode 
    switch bapGeometry                                                                                        % switch/case to assign insignificant regions as no-flux
        case 1                                                                                                % case 1 - planar slab
            noFluxEnds.selection.named('geom1_csel2_bnd');                                                    % assign the hydrogel slab side faces to the no flux subnode
        case 2                                                                                                % case 2 - cylidner
            noFluxEnds.selection.named('geom1_csel3_bnd');                                                    % assign hydrogel cylinder end faces to the no-flux subnode
        case 3                                                                                                % case 3 - hollow cylinder
            noFluxEnds.selection.named('geom1_csel3_bnd');                                                    % assign hydrogel cylinder end faces to the no-flux subnode
    end                                                                                                       % end switch/case


    %% 5.10. Create stationary study and solution node
    % 5.10.1. Create a study
    model.study.create('std1');
    model.study('std1').create('stat', 'Stationary');
    % 5.10.2. Create a solution
    model.sol.create('sol1');
    model.sol('sol1').study('std1');
    model.sol('sol1').attach('std1');
    model.sol('sol1').create('st1', 'StudyStep');
    model.sol('sol1').create('v1', 'Variables');
    model.sol('sol1').create('s1', 'Stationary');
    model.sol('sol1').feature('s1').create('fc1', 'FullyCoupled');
    model.sol('sol1').feature('s1').create('i1', 'Iterative');
    model.sol('sol1').feature('s1').feature('i1').create('mg1', 'Multigrid');
    model.sol('sol1').feature('s1').feature('i1').feature('mg1').feature('pr').create('sl1', 'SORLine');
    model.sol('sol1').feature('s1').feature('i1').feature('mg1').feature('po').create('sl1', 'SORLine');
    model.sol('sol1').feature('s1').feature('i1').feature('mg1').feature('cs').create('d1', 'Direct');
    model.sol('sol1').feature('s1').feature.remove('fcDef');
    % 5.10.3. Create plot groups and numerical data sets 
    model.result.dataset.create('av1', 'Average');
    model.result.dataset('av1').selection.named('geom1_csel1_dom');
    model.result.numerical.create('pev1', 'EvalPoint');
    model.result.numerical('pev1').set('probetag', 'none');
    model.result.create('pg1', 'PlotGroup3D');
    model.result.create('pg2', 'PlotGroup3D');
    model.result('pg1').create('slc1', 'Slice');
    model.result('pg2').create('surf1', 'Surface');
    % 5.10.4. Attach solution to study 1 and set defaults
    model.sol('sol1').attach('std1');
    model.sol('sol1').feature('s1').feature('fc1').set('initstep', 0.01);
    model.sol('sol1').feature('s1').feature('fc1').set('minstep', 1.0E-6);
    model.sol('sol1').feature('s1').feature('fc1').set('maxiter', 50);
    model.sol('sol1').feature('s1').feature('i1').set('nlinnormuse', true);
    model.sol('sol1').feature('s1').feature('i1').set('maxlinit', 400);
    model.sol('sol1').feature('s1').feature('i1').set('rhob', 40);
    model.sol('sol1').feature('s1').feature('i1').feature('mg1').feature('pr').feature('sl1').set('linerelax', 0.2);
    model.sol('sol1').feature('s1').feature('i1').feature('mg1').feature('pr').feature('sl1').set('relax', 0.4);
    model.sol('sol1').feature('s1').feature('i1').feature('mg1').feature('po').feature('sl1').set('linerelax', 0.2);
    model.sol('sol1').feature('s1').feature('i1').feature('mg1').feature('po').feature('sl1').set('seconditer', 2);
    model.sol('sol1').feature('s1').feature('i1').feature('mg1').feature('po').feature('sl1').set('relax', 0.4);
    model.sol('sol1').feature('s1').feature('i1').feature('mg1').feature('cs').feature('d1').set('linsolver', 'pardiso');
    model.sol('sol1').feature('s1').feature('i1').feature('mg1').feature('cs').feature('d1').set('pivotperturb', 1.0E-13);

    
    %% 5.11. Run study
    try                                                                                   % try/catch statement to catch solving errors    
        model.sol('sol1').runAll;                                                         % run the study   
    catch                                                                                 % If there's an error...
        didFail = 1;                                                                      % ...update the failure indicator
        solverFailureMessage = ['\nError: Failed to solve. Restarting iteration now.\n']; % write solver failure message
        fprintf(solverFailureMessage)                                                     % print failure message
        clearvars xPos yPos zPos isletDiameters isletVolumes                              % clear islet variables
        continue                                                                          % ... return to the top of the while loop   
    end                                                                                   % end try/catch statement

    
    %% 5.12. Collect population data
    % 5.12.1. Create export tables
    model.result.table.create('tbl1', 'Table');                           % make table 1 for popAvpO2
    model.result.table.create('tbl2', 'Table');                           % make table 2 for NVF
    model.result.table.create('tbl3', 'Table');                           % make table 3 for S
    % 5.12.2. Average volume in islets 
    model.result.numerical.create('av1', 'AvVolume');                     % create average volume node 
    model.result.numerical('av1').selection.named('geom1_csel1_dom');     % select islets
    model.result.numerical('av1').set('probetag', 'none');                % name the probe nothing
    model.result.numerical('av1').set('table', 'tbl1');                   % set result to table 1
    model.result.numerical('av1').set('expr', {'pO2'});                   % define calculation (average pO2)
    model.result.numerical('av1').set('unit', {'mmHg'});                  % define unit 
    model.result.numerical('av1').set('descr', {''});                     % add blank description
    model.result.numerical('av1').setResult;                              % calculate value, and send result to table 1
    % 5.12.3. Necrotic volume
    model.result.numerical.create('int1', 'IntVolume');                   % create volume integration node 
    model.result.numerical('int1').selection.named('geom1_csel1_dom');    % select islets
    model.result.numerical('int1').set('probetag', 'none');               % name the probe nothing
    model.result.numerical('int1').set('table', 'tbl2');                  % set result to table 2
    model.result.numerical('int1').set('expr', {'pO2 < (0.081 [mmHg])'}); % define threshold representing necrotic tissue
    model.result.numerical('int1').set('unit', {'m^3'});                  % define unit
    model.result.numerical('int1').set('descr', {''});                    % add blank description
    model.result.numerical('int1').setResult;                             % calculate value, and send result to table 2
    % 5.12.4. Insulin secretion potential
    model.result.numerical.create('av2', 'AvVolume');                     % create volume integration node 
    model.result.numerical('av2').selection.named('geom1_csel1_dom');     % select islets 
    model.result.numerical('av2').set('probetag', 'none');                % name the probe nothing
    model.result.numerical('av2').set('table', 'tbl3');                   % set result to table 3
    model.result.numerical('av2').set('expr', {'1-sIns'});                % define threshold representing anoxic tissue
    model.result.numerical('av2').set('unit', {''});                      % define unit (unitless)
    model.result.numerical('av2').set('descr', {''});                     % add blank description
    model.result.numerical('av2').setResult;                              % calculate value, and send result to table 3
    % 5.12.5. Collect Data in MATLAB arrays
    popAvpO2(ind) = model.result.table('tbl1').getReal;                   % collect popAvpO2 data in matlab array
    popNVF(ind) = model.result.table('tbl2').getReal/sum(isletVolumes);   % collect popNVF data in matlab array
    popPsi(ind) = model.result.table('tbl3').getReal;                     % collect popS data in matlab array

    
    %% 5.13. Save model
    switch ind                                                           % switch/case to only save the first iteration's model
        case 1                                                           % if it is the first model
            modelFileName = [sharpFolder '\COMSOL model\YourModel.mph']; % write model file name
            mphsave(model, modelFileName)                                % save model file
    end                                                                  % end switch/case
    
    %% 5.14. Calculate inidividual iselt outcomes
    % 5.14.1. Mean pO2
    i1 = 1;                                                                                   % indexer
    isletAvpO2 = zeros(size(isletDiameters));                                                 % preallocate isletAvpO2 array   
    for i1 = 1:length(isletDiameters)                                                         % loop over all islets
        avTag1 = model.result.numerical.uniquetag('av');                                      % create a unique average tag
        isletAvpO2VolNode(i1) = model.result.numerical.create(avTag1, 'AvVolume');            % create a unique average volume node
        selectionname = ['geom1_' 'sph' num2str(i1) '_dom'];                                  % update selection name with islet number
        isletAvpO2VolNode(i1).selection.named(selectionname);                                 % set the unique average volume node to target indicated sphere
        isletAvpO2VolNode(i1).set('probetag', 'none');                                        % set no probetag
        tbltag1 = model.result.table.uniquetag('tbl');                                        % create a unique table tag
        model.result.table.create(tbltag1, 'Table');                                          % create a unique table
        isletAvpO2VolNode(i1).set('table', tbltag1);                                          % set the table for the average volume node to the unique table
        isletAvpO2VolNode(i1).set('expr', {'pO2'});                                           % define calculation (average pO2)
        isletAvpO2VolNode(i1).set('unit', {'mmHg'});                                          % define unit                
        isletAvpO2VolNode(i1).set('descr', {''});                                             % add blank description
        isletAvpO2VolNode(i1).setResult;                                                      % calculate value, and send it to unique table
        isletAvpO2(i1) = model.result.table(tbltag1).getReal;                                 % collect isletAvpO2 data in matlab array
    end                                                                                       % end islet data extraction loop
    cellIsletAvpO2{ind} = isletAvpO2;                                                         % collect islet pO2 for particular iteration in a cell
    % 5.14.2. Loss of insulin secretion capacity (Psi)
    i2 = 1;                                                                                   % indexer
    isletPsi = zeros(size(isletDiameters));                                                   % preallocate loss of function array   
    for i2 = 1:length(isletDiameters)                                                         % loop over all islets
        avTag2 = model.result.numerical.uniquetag('av');                                      % create a unique average tag
        isletAVFNode(i2) = model.result.numerical.create(avTag2, 'AvVolume');                 % create a unique average volume node
        selectionname = ['geom1_' 'sph' num2str(i2) '_dom'];                                  % update selection name with islet number
        isletAVFNode(i2).selection.named(selectionname);                                      % set the unique average volume node to target indicated sphere
        isletAVFNode(i2).set('probetag', 'none');                                             % set no probetag
        tbltag2 = model.result.table.uniquetag('tbl');                                        % create a unique table tag
        model.result.table.create(tbltag2, 'Table');                                          % create a unique table
        isletAVFNode(i2).set('table', tbltag2);                                               % set the table for the average volume node to the unique table
        isletAVFNode(i2).set('expr', {'1-sIns'});                                             % define calculation (Psi)
        isletAVFNode(i2).set('unit', {''});                                                   % define unit                
        isletAVFNode(i2).set('descr', {''});                                                  % add blank description
        isletAVFNode(i2).setResult;                                                           % calculate value, and send it to unique table
        isletPsi(i2) = model.result.table(tbltag2).getReal;                                   % collect isletAVF data in matlab array
    end                                                                                       % end islet data extraction loop
    cellIsletS{ind} = isletPsi;                                                               % collect islet Psi for particular iteration in a cell
    % 5.14.3. Necrosis
    i3 = 1;                                                                                   % indexer
    isletNVF = zeros(size(isletDiameters));                                                   % preallocate isletAvpO2 array   
    for i3 = 1:length(isletDiameters)                                                         % loop over all islets
        avTag3 = model.result.numerical.uniquetag('av');                                      % create a unique average tag
        isletNVFNode(i3) = model.result.numerical.create(avTag3, 'AvVolume');                 % create a unique average volume node
        selectionname = ['geom1_' 'sph' num2str(i3) '_dom'];                                  % update selection name with islet number
        isletNVFNode(i3).selection.named(selectionname);                                      % set the unique average volume node to target indicated sphere
        isletNVFNode(i3).set('probetag', 'none');                                             % set no probetag
        tbltag3 = model.result.table.uniquetag('tbl');                                        % create a unique table tag
        model.result.table.create(tbltag3, 'Table');                                          % create a unique table
        isletNVFNode(i3).set('table', tbltag3);                                               % set the table for the average volume node to the unique table
        isletNVFNode(i3).set('expr', {'pO2 < (0.081 [mmHg])'});                               % define calculation (NVF)
        isletNVFNode(i3).set('unit', {'m^3'});                                                % define unit                
        isletNVFNode(i3).set('descr', {''});                                                  % add blank description
        isletNVFNode(i3).setResult;                                                           % calculate value, and send it to unique table
        isletNVF(i3) = model.result.table(tbltag3).getReal;                                   % collect isletNVF data in matlab array
    end                                                                                       % end islet data extraction loop
    cellIsletNVF{ind} = isletNVF;                                                             % collect islet NVF for particualr iteration in a cell
    % 5.14.4. Collate data
    cellIsletInfo{ind} = [xPos' yPos' zPos' isletDiameters' isletAvpO2' isletNVF' isletPsi']; % collated x, y, z poswitions, islet diameter, and outcomes for all islets for each iteration
    saveFile = cellIsletInfo{ind};                                                            % isolate individual iteration results
    
    
    %% 5.15. Save data to file
    individualIsletData = [sharpFolder '\Individual islet data\isletDataIteration_' num2str(ind) '.txt']; % write file name, attaching iteration indexer at the end
    save(individualIsletData, 'saveFile')                                                                 % save the file

    
    %% 5.16. Clear COMSOL model and reset failure indicator
    clearvars model geom1 % clear model and geometry
    didFail = 0;          % reset failure indicator
    
    
end % end study loop 


%% 6. Print islet population results
netResults = table(popAvpO2',popPsi',popNVF') % column 1: avpO2; column 2: loss of insulin secretion capacity, column 3: necrosis


%% 7. Audibly indicate completion
load handel % load sound
sound(y,Fs) % make sound
