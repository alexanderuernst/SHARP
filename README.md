% Author: Alexander Ulrich Ernst, PhD (email: aue2@cornell.edu)
% Date created: 20Jul20
% Last updated: 28Jul22

% This program calculates oxygen distributions and related outcomes in user
% -programmable bioartificial pancreas devices comprising hydrogel and 
% islets.

% INSTRUCTIONS
% Follow input prompts to enter device contents: cell type, cell volume, 
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

% If the user wishes to change the values of individual parameters (e.g., 
% the value of the external oxygen concentrioton), this can be done 
% manually in section 5.6. wherein the parameter values are hard-coded. The 
% default path is set to the author's Dropbox, but can (and should) be 
% changed manually by the user.
