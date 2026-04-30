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
% load '\\david.lse.thm.de\Jessica\MATLAB\phantoms\Patient_0_Flatten.mat'
load '\\david.lse.thm.de\Jessica\MATLAB\phantoms\Waterphantom_10x10x32cm_target_5cm.mat'
% % Show an image of the CT slice with the lung:
% image(ct.cubeHU{1,1}(:,:,50),'CDataMapping','scaled');
% xlabel('x in mm');
% ylabel('y in mm');
% colorbar


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

% Test worst case: Pmod lokal = Pmod konstant!
% ct.cube_pmod(ct.cube_pmod~=0)=750;
%% 

% Definition of the block size
block_size = [cluster_size, cluster_size, cluster_size];


pln.radiationMode   = 'carbon';     % either photons / protons / carbon
pln.machine         = 'HIT_APM';

modelName           = 'none';
quantityOpt         = 'physicalDose';   

% The remaining plan parameters are set like in the previous example files
pln.numOfFractions = 1;

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


%% Recalculation with MIT basedata
pln.machine         = 'MIT_carbon';

% Basedata HIT
load('C:\Users\Stolzenberg\Documents\MATLAB\matRad_official\basedata\carbon_HIT_ORG.mat')
machine_HIT = machine;

% Basedata MIT
load('\\david.lse.thm.de\Jessica\MATLAB\basedata\carbon_MIT_carbon.mat')

% MIT-Energien
MIT_energies = [machine.data.energy];
HIT_energies = [machine_HIT.data.energy];

% Schleife über alle Rays
stf_MIT = stf;
for r = 1:length(stf.ray)
    energies = stf.ray(r).energy;           % Energien dieser Ray
    nearest_val = zeros(size(energies));    % nächster Wert in MIT

    for i = 1:length(energies)
        [~, idx] = min(abs(MIT_energies - energies(i)));
        nearest_val(i) = MIT_energies(idx);
    end

    % Ergebnisse speichern
    stf_MIT.ray(r).energy = nearest_val;
end
stf_MIT.SAD = machine.meta.SAD;
%% hier die Koordinaten der Rays und der stf.sourcePoint und stf.sourcePoint_bev ändern!
stf_MIT.sourcePoint_bev = [0 -stf_MIT.SAD 0];


% get (active) rotation matrix
% transpose matrix because we are working with row vectors
rotMat_vectors_T = transpose(matRad_getRotationMatrix(pln.propStf.gantryAngles,pln.propStf.couchAngles));

stf_MIT.sourcePoint = stf_MIT.sourcePoint_bev*rotMat_vectors_T;

% for j=1:length(stf_MIT.ray)
%     stf_MIT.ray(j).targetPoint_bev = [stf_MIT.ray(j).targetPoint_bev(1) ...
%         stf_MIT.SAD ...
%         stf_MIT.ray(j).targetPoint_bev(3)];
%     stf_MIT.ray(j).targetPoint = stf_MIT.ray(j).targetPoint_bev*rotMat_vectors_T;
%     %% Bestimmung des focusIx:
%     currentMinimumFWHM = matRad_interp1(machine.meta.LUT_bxWidthminFWHM(1,:)',...
%         machine.meta.LUT_bxWidthminFWHM(2,:)',...
%         pln.propStf.bixelWidth, ...
%         machine.meta.LUT_bxWidthminFWHM(2,end));
%     focusIx  =  ones(stf_MIT.numOfBixelsPerRay(j),1);
%     [~, vEnergyIx] = min(abs(bsxfun(@minus,[machine.data.energy]',...
%         repmat(stf_MIT.ray(j).energy,length([machine.data]),1))));
% 
%     % get for each spot the focus index
%     for k = 1:stf_MIT.numOfBixelsPerRay(j)
%         focusIx(k) = find(machine.data(vEnergyIx(k)).initFocus.SisFWHMAtIso > currentMinimumFWHM,1,'first');
%     end
% 
%     stf_MIT.ray(j).focusIx = focusIx';
% end 

% %% neue Herangehensweise für die Änderung des Focus:
% % Finde die SisFWHM, die an die Focusstufe vom HIT am nächsten kommt:
% for j=1:length(stf_MIT.ray)
%     focusIx = zeros([length(stf_MIT.ray(j).energy),1]);
%     for k=1:length(stf_MIT.ray(1).energy)
%         [~, idx_EMIT] = min(abs(MIT_energies - stf_MIT.ray(j).energy(k)));
%         [~, idx_EHIT] = min(abs(HIT_energies - stf.ray(j).energy(k)));
%         % Finde SisFWHM passend zur Fokusstufe in HIT Basisdaten
%         [~, idx] = min(abs(machine.data(idx_EMIT).initFocus.SisFWHMAtIso - machine_HIT.data(idx_EHIT).initFocus.SisFWHMAtIso(stf.ray(j).focusIx(k))));
%         focusIx(k) = idx;
%     end
%     % Ergebnisse speichern
%     stf_MIT.ray(j).focusIx = focusIx';
% end

%% Anpassung des Focus an die klinisch benutzbaren Werte vom MIT:
for j=1:length(stf_MIT.ray)
    focusIx = zeros([length(stf_MIT.ray(j).energy),1]);
    for k=1:length(stf_MIT.ray(1).energy)
        % Finde SisFWHM passend zur klinischen Fokusstufe in MIT Basisdaten
        if stf_MIT.ray(j).energy(k) >= 224.22
            focusIx(k) = 3;
        elseif stf_MIT.ray(j).energy(k) >= 164.09 && stf_MIT.ray(j).energy(k) <= 222.97
            focusIx(k) = 2;
        elseif stf_MIT.ray(j).energy(k) <=162.56
            focusIx(k) = 1;
        end
        
    end
    % Ergebnisse speichern
    stf_MIT.ray(j).focusIx = focusIx';
end
%% Anpassen der weights an die unterschiedlichen Peakdosen zwischen 
% HIT und MIT in den Basisdaten:
% n=1;
% for j = 1 : dij.numOfRaysPerBeam
%     for k = 1 : length(stf.ray(1).energy)
%         [~, idx_EMIT] = min(abs(MIT_energies - stf_MIT.ray(j).energy(k)));
%         [~, idx_EHIT] = min(abs(HIT_energies - stf.ray(j).energy(k)));
%         % Finde SisFWHM passend zur Fokusstufe in HIT Basisdaten
%         wcorr = max(machine_HIT.data(idx_EHIT).Z.doseORG)/max(machine.data(idx_EMIT).Z);
%         resultGUI_homogeneous_MIT.w(n) = resultGUI_homogeneous.w(n)*wcorr;
%         n=n+1;
%     end
% end

% % Test mit stf Generierung der eigentlichen Basisdaten
% stf_MIT2 = matRad_generateStf(ct,cst,pln);

% Nachrechnen des HIT Plans mit MIT Basisdaten:
resultGUI_homogeneous_MIT = matRad_calcDoseDirect(ct,stf_MIT,pln,cst,resultGUI_homogeneous.w);

%% Compare dose of MIT and HIT plan:
% filedir = "\\david.lse.thm.de\Jessica\MATLAB\Results\Vergleich_MIT_HIT\C12";
% lungwing_tumor = "Lung L";
% compareDose_plots(resultGUI_homogeneous.physicalDose,resultGUI_homogeneous_MIT.physicalDose,ct,cst,[1 1 0], [], pln, [1,1], 0,'local', filedir, "HIT__vs_MIT" );
% quality_indices = DVH_analysis(resultGUI_homogeneous.physicalDose,resultGUI_homogeneous_MIT.physicalDose, ct, cst, pln, filedir, "HIT", "MIT", lungwing_tumor);

ct.planDose = 60;
[gammaCube,gammaPassRate,hfig] = matRad_compareDose(resultGUI_homogeneous.physicalDose./max(resultGUI_homogeneous.physicalDose(:)),resultGUI_homogeneous_MIT.physicalDose./max(resultGUI_homogeneous_MIT.physicalDose(:)), ct, cst,[1 1 1], 'off',pln, [3,3],0, 'global');
% [gammaCube,gammaPassRate,hfig] = matRad_compareDose(resultGUI_homogeneous.physicalDose,resultGUI_homogeneous_MIT2.physicalDose, ct, cst,[1 1 1], 'off',pln, [3,3],0, 'global');

%%
% cst_withLungFlag = pln.propHeterogeneity.cstHeteroAutoassign(cst);

% %% Calculate dose again with heterogeneityCorrection
% % Local:
% pln.propHeterogeneity = matRad_HeterogeneityConfig();
% pln.propHeterogeneity.calcHetero = 1;
% resultGUI_heterogeneous_Flatten = matRad_calcDoseDirect(ct,stf,pln,cst_withLungFlag,resultGUI_homogeneous.w);
% 
% % Constant Pmod:
% pln.propHeterogeneity.type = 'voxelwise';
% pln.propHeterogeneity.modPower = mean(pmodCT.cube_pmod(pmodCT.cube_pmod(:)~=0)); %% Mittelwert der Pmod aus dem lokalen PmodCT!
% % pln.propHeterogeneity.modPower = 750;
% pln.propHeterogeneity.calcHetero = 1;
% 
% resultGUI_heterogeneous_Winter = matRad_calcDoseDirect(ct,stf,pln,cst_withLungFlag,resultGUI_homogeneous.w);

%% Comparison of Models:
% matRad_compareDose(carbHomo.physicalDose,carbHetero.physicalDose,ct,cst,[1 0 0]);
% Homo vs. Flatten
% matRad_compareDose(resultGUI_homogeneous.physicalDose,resultGUI_heterogeneous_Flatten.physicalDose,ct,cst,[1 1 0]);
% Flatten vs. Winter
% matRad_compareDose(resultGUI_heterogeneous_Flatten.physicalDose,resultGUI_heterogeneous_Winter.physicalDose,ct,cst,[1 1 0], [], pln, [3,3], 0,'local');

