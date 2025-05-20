% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% Created by: Jessica Stolzenberg on 16.01.2025, last updated: 20.01.2025
% Purpose: Implement model by Flatten et al. (2021) to a CT
% Input: CT with contours of the lung, each lung side is handled
% individually
% These parameters can be changed manually for the cluster process:
% cluster_size (size of the resulting blocks) and edgemethod
% edgemethod: 'None', 'Simple' to include the edges for the calculation of
% Pmod
% Steps: 
% 0) Load CT
% 1) Extract the lungs out of the original CT (CT with HU values)
% 2) Cluster the voxels
% 3) Perform a gaussian fit
% 4) Insert the sigma and mu values into the function by Flatten et
% al.(2021)
% 5) Set all values of each Cluster to the same Pmod value.
% 6) Save the Pmod valuemap and set all other values to zero
% Output: CT with modulation powers instead of the HU or density values in
% the lung, everywhere else Pmod should be 0.
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
tic
clear; clc
%% 0) Load CT
load '\\david.lse.thm.de\Jessica\MATLAB\phantoms\Patient_0_Flatten.mat'

% Show an image of the CT slice with the lung:
image(ct.cubeHU{1,1}(:,:,50),'CDataMapping','scaled');
xlabel('x in mm');
ylabel('y in mm');
colorbar


cluster_size = 7; % Changes the size of the resulting blocks
bin_division = 10;
order = 'zyx';
edgemethod = 'Local';
%% Test:
% load('\\david.lse.thm.de\Jessica\MATLAB\Cluster_Algorithm_Marc\Cluster_Robustness.mat')
% pmodCT = cs_ct(2);
% ct.cube_pmod = pmodCT.pmod;

load('C:\Users\Stolzenberg\Documents\MATLAB\matRad_official\Stolzenberg_files\pmodCT_cs_7_ohne_HU_range.mat')
% [pmodCT] = pmodBloxx(ct,cst, cluster_size, bin_division, order, edgemethod);
ct.cube_pmod = pmodCT.cube_pmod;

% Test worst case:
ct.cube_pmod(ct.cube_pmod~=0)=750;
%% Test:
% ct.cube_pmod = 250;


% Definition of the block size
block_size = [cluster_size, cluster_size, cluster_size];


pln.radiationMode   = 'protons';     % either photons / protons / carbon
pln.machine         = 'Generic_APM';

modelName           = 'none';
quantityOpt         = 'physicalDose';   

% The remaining plan parameters are set like in the previous example files
pln.numOfFractions = 6;

pln.propStf.gantryAngles  = 0;
pln.propStf.couchAngles   = 0;
pln.propStf.bixelWidth    = 10;
pln.propStf.numOfBeams    = numel(pln.propStf.gantryAngles);
pln.propStf.isoCenter     = ones(pln.propStf.numOfBeams,1) * matRad_getIsoCenter(cst,ct,0);
pln.propOpt.runDAO        = 0;
pln.propOpt.runSequencing = 0;
%% Modification: DoseGridresolution = ct resolution
% pln.propDoseCalc.doseGrid.resolution = ct.resolution;
%%
% retrieve bio model parameters
pln.bioParam = matRad_BioModel(pln.radiationMode,quantityOpt,modelName);

% retrieve scenarios for dose calculation and optimziation
pln.multScen = matRad_multScen(ct,'nomScen'); % optimize on the nominal scenario                                            

pln.propHeterogeneity = matRad_HeterogeneityConfig();
pln.propHeterogeneity.calcHetero = 1;
%% Generate Beam Geometry STF
stf = matRad_generateStf(ct,cst,pln);
% stf = matRad_generateStfPencilBeam(pln,ct);

%%
dij = matRad_calcParticleDose(ct,stf,pln,cst);

%% Inverse Optimization  for IMPT based on RBE-weighted dose
resultGUI_homogeneous = matRad_fluenceOptimization(dij,cst,pln);

cst_withLungFlag = pln.propHeterogeneity.cstHeteroAutoassign(cst);

%% Calculate dose again with heterogeneityCorrection
pln.propHeterogeneity = matRad_HeterogeneityConfig();
pln.propHeterogeneity.calcHetero = 1;
resultGUI_heterogeneous_Flatten = matRad_calcDoseDirect(ct,stf,pln,cst_withLungFlag,resultGUI_homogeneous.w);

pln.propHeterogeneity.type = 'voxelwise';
pln.propHeterogeneity.modPower = 750;
pln.propHeterogeneity.calcHetero = 1;

resultGUI_heterogeneous_Winter = matRad_calcDoseDirect(ct,stf,pln,cst_withLungFlag,resultGUI_homogeneous.w);

%% Visualize differences
% matRad_compareDose(carbHomo.physicalDose,carbHetero.physicalDose,ct,cst,[1 0 0]);
% Homo vs. Flatten
% matRad_compareDose(resultGUI_homogeneous.physicalDose,resultGUI_heterogeneous_Flatten.physicalDose,ct,cst,[1 1 0]);
% Flatten vs. Winter
matRad_compareDose(resultGUI_heterogeneous_Flatten.physicalDose,resultGUI_heterogeneous_Winter.physicalDose,ct,cst,[1 1 0], [], pln);


%  Plot parameter:
set(groot, 'defaultLineLineWidth', 3);          % Standard-Linienstärke
set(groot, 'defaultAxesFontSize', 16);         % Standard-Schriftgröße der Achsen
set(groot, 'defaultAxesLineWidth', 1.5);       % Rahmenbreite der Achsen
set(groot, 'defaultTextFontSize', 18);         % Schriftgröße für Titel und Beschriftungen
set(groot, 'defaultFigureColor', 'w');         % Hintergrundfarbe der Figur (weiß)
set(groot, 'defaultAxesGridLineStyle', '--');  % Stil der Gitternetzlinien
set(groot, 'DefaultTextInterpreter', 'latex');
set(groot, 'DefaultAxesTickLabelInterpreter', 'latex');
set(groot, 'DefaultLegendInterpreter', 'latex');
set(groot, 'DefaultColorbarTickLabelInterpreter', 'latex'); 
toc