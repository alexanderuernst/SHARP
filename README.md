**Simulated Heterogenity and Randomness Program (SHARP)**

Alexander Ulrich Ernst, PhD (email: aue2@cornell.edu)

September 2, 2022

**1. Description**

This program calculates oxygen distributions and related outcomes in user-programmable bioartificial pancreas devices comprising hydrogel and islets. Data collected using this program is described in an article entitled _A predictive computational platform for optimizing the design of biartificial pancreas devices_, currently under review at _Nature Communications_.


**2. Program Requirements and Installation**

This program was developed in COMSOL Multiphyusics 5.4 with MATLAB (v. 2019a). It should be usable with any version of COMSOL Multiphysics and MATLAB, though the "COMSOL Multiphysics with MATLAB" portal, which can be purchased from COMSOL is required. To install the program, download the "SHARP Application" folder. Then, launch the COMSOL Multiphysics with MATLAB portal, open the "sharApplicationv1" script, and set the path to the location of the downloaded folder. Finally, copy the path in the first and only line of _Section 3_. The program was developed and exclusively operated on a Lenovo X1 Carbon with a RAM of 15.8 GB and a clockspeed of 2.11 GHz running Mcirosoft WIndows 11 operating system. 


**3. Insructions**

Follow input prompts to enter device contents: cell type, cell volume, device type, device parameters, and number of Monte Carlo iterations you wish the program to run. A minimum of three iterations, or a total islet volume of 1500 IEQ, is recommended. The program will calculate three device fates: the mean pO2, net loss of insulin secretory capacity (Psi), and net volumetric fraction of necrotic tissue of the islets. It will calculate these outcomes on a population basis (e.g., the average pO2 of all the encapsulated islets) and on an individual islet basis. The former is printed as a table and the latter is saved as matlab files to a designated folder. In addition, the program will display the device geometry for each iteration of the simulation and will save the first iteration as a COMSOL file for the user's convenience. For each iteration of the simulation, a readout of the device specifications is be provided for user confirmaton. 


**4. Instructions for Custom Use**

If the user wishes to change the values of individual parameters (e.g., the value of the external oxygen concentrioton), this can be done manually in section 5.6. wherein the parameter values are hard-coded. If the user wishes to change values of islet properties, namely the size distribution parameters and the oxygen consumption rate, the user may modify constants in _Section 2.2._ of the function file "sharpSeedingAlgo".


**5. Example Test Case**

To confirm that the program is functioning properly, I recommend that the following inputs are provided (estimated run time: 5 minutes and 40 seconds):


  _Enter BAP geometry type (1 = planar slab, 2 = cylinder, 3 = hollow cylinder):_  1
  
  _Enter islet source (1 = human, 2 = NN SCBs, 3 = MM SCBs, 4 = pig, 5 = rat, 6 = mouse):_  2
  
  _Enter the total islet volume, in IEQ, encapsulated in the BAP device:_  50
  
  _Enter the thickness of the hydrogel slab (in micrometers):_  500
  
  _Enter the diameter of the hydrogel slab (in millimeters):_  5
  
  _Enter number of iterations you wish the Monte Carlo simulation to perform:_  3

Values resembling the following output should be shown in the Command Window of the MATLAB workspace (note, the values will not be exactly the same as the output of the results are stochastic):

![image](https://user-images.githubusercontent.com/91346875/188233828-0362ab09-11c5-4006-b2f4-4f1fddd9a4e0.png)

Above, Var1, Var2, and Var3 represent the mean pO2, net loss of function, and net necrosis of the islet population, respectively.

In addition, a COMSOL file entitled "YourModel.m" should appear in the "COMSOL Model" subfolder, and three MATLAB files entitled "isletDataIteration_1", "isletDataIteration_2" and "isletDataIteration_3" should appear in the "Individual islet data" subfolder.


**6. Reproducing Manuscript Results**

Results near those provided in the manuscript _A predictive computational platform for optimizing the design of bioartificial pancreas devices_, currently under review at _Nature Communications_, can be reproduced by inputting the parameters of the devices as described in the manuscript text and supplementary information. Again, exact values cannot be obtained as the program results are stochastic. 
