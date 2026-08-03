% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% Created by: Jessica Stolzenberg on 07.11.2025, last updated: 07.11.2025
% Purpose: Patient study with five patients using pmodBloxx
% Input: CT with contours of the lung, each lung side is handled
% individually
% Steps:
% 0) Load CT
% 1) Calculate homogenous plan
% 2) Calculate heterogeneous plan using pmodBloxx
% 3) Compare both plans
% Output: resultGUI homogeneous and heterogeneous calculated using
% pmodBloxx
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
tic
clear; clc
%% CT file paths:
 Patients = ["Patient_0_Flatten.mat"];%, "Patient_2_Flatten.mat", "Patient_4_Flatten.mat"];
% Patients = ["Annika1.mat"];
folder = 'C:\Users\Stolzenberg\Documents\MATLAB\matRad_official\phantoms';
folder_Pmod = '\\david.lse.thm.de\Jessica\MATLAB\Pmod_Patients';
filedir_gamma = fullfile(folder_Pmod, "Gamma");
filedir=fullfile(folder_Pmod, "DVH");
filedir_pvh=fullfile(folder_Pmod, "PVH");


step_size = 0.01; % Parameter for DVH    
step_size_pmod = 10; % Parameter for DVH    

%% default options for pmodBloxx
cluster_size = 10; % Changes the size of the resulting blocks
bin_division = 12;
order = 'zyx';
edgemethod = 'Local';


%% Loop over all patients:
for i = 1:length(Patients)
    path = fullfile(folder, Patients(i));
    load(char(path))
    % token = regexp(Patients(i), '^(Patient_\d+)', 'tokens');
    % patientID = token{1}{1};
    patientID = "2";
    patientTitle = strrep(patientID, '_', ' ');
    % cst{17,6}=[];
    %% dose constraints des PTV wie in Marburg (nach Kilian) --> D_min = 100%, D_max = 100% und D_max = 105%
    % hier werden für die Pläne bisher nur Dose constraints für das PTV
    % verwendet, da andere Constraints zu schlechten Plänen geführt haben,
    % d.h. im DVH konnte man sehen, dass das PTV nicht 100% mit Dosis
    % versorgt wurde, sowie starke Abweichungen in den Idealwerten des
    % Konformitätsindex und Homogenitätsindex --> kann man ändern!
    ptv = {'PTV'}; ctv = {'CTV'}; body = {'BODY'};
    try
        ptvIdx = find(contains(cst(:,2), ptv));
        ctvIdx = find(contains(cst(:,2), ctv));
        bodyIdx = find(contains(cst(:,2), body));
        cst{ctvIdx,6} = [];
        cst{ptvIdx,6} = [];
        cst{ptvIdx,6}{1} = struct(DoseObjectives.matRad_SquaredUnderdosing(40,30)); % Dmin=100%, 40
        % cst{ptvIdx,6}{2} = struct(DoseObjectives.matRad_SquaredOverdosing(1000,30)); % Dmax = 100%, 4 (raus?) nur PTV
        cst{ptvIdx,6}{2} = struct(DoseObjectives.matRad_SquaredOverdosing(40,30.5)); % Dmax = 105%, 40
        % cst{ctvIdx,6}{1} = struct(DoseObjectives.matRad_SquaredUnderdosing(1000,30)); % Dmin=100%
        % cst{ctvIdx,6}{1} = struct(DoseObjectives.matRad_SquaredOverdosing(1000,31.5)); % Dmax = 100% 
        % cst{ctvIdx,6}{3} = struct(DoseObjectives.matRad_SquaredOverdosing(1000,31.5)); % Dmax = 105% 
        % cst{bodyIdx,6}{1} = struct(DoseObjectives.matRad_SquaredOverdosing(100,5)); % Dmax = 5 Gy raus
    catch
        warning('PTV/CTV not found for Conformity Index calculation.');
    end

    %% calculate PmodCT for patient:
    % Load PmodCT files for patient:
    pmod_path = fullfile(folder_Pmod, append('Pmod_',Patients(i)));
    load(char(pmod_path))
    pmodct.pmod(pmodct.pmod<0) = 0;
    pmodct.pmod(pmodct.pmod>900) = 900;
    pmodCT(i).pmod = pmodct.pmod;
    
    %% Den Teil auskommentieren für Lungenmodulation weglassen:
    % Calculate pmodCT locally inside matLab:
    % [pmodCT_patient] = pmodBloxx(ct,cst, cluster_size, bin_division, order, edgemethod);
    % Halte Pmod in range von 0um bis 900um:
    
    % pmodCT(i).pmod = pmodCT_patient.cubepmod;
    %%
    % Definition of the block size
    block_size = [cluster_size, cluster_size, cluster_size];


    pln.radiationMode   = 'carbon';     % protons / carbon
    pln.machine         =  'HIT_APM';% 'Generic_APM';

    %modelName           = 'none';
    %quantityOpt         = 'physicalDose';

    % modelName    = 'constRBE';
    % quantityOpt  = 'RBExD';

    modelName    = 'LEM';
    quantityOpt  = 'RBExD';


    % The remaining plan parameters are set like in the previous example files
    pln.numOfFractions = 10;

    pln.propStf.gantryAngles  = 90;
    pln.propStf.couchAngles   = 0;
    pln.propStf.bixelWidth    = 2; % 1.8 mm für C12, 2/3 mm für Protonen (MIT)
    pln.propStf.numOfBeams    = numel(pln.propStf.gantryAngles);
    pln.propStf.isoCenter     = ones(pln.propStf.numOfBeams,1) * matRad_getIsoCenter(cst,ct,0);
    pln.propOpt.runDAO        = 0;
    pln.propOpt.runSequencing = 0;
    %% Modification: DoseGridresolution = ct resolution
    % pln.propDoseCalc.doseGrid.resolution = ct.resolution;
    % Default: Dose Grid Resolution = [3 3 3] (wie in Studie aus Japan)
    %%
    % retrieve bio model parameters
    pln.bioParam = matRad_BioModel(pln.radiationMode,quantityOpt,modelName);

    % retrieve scenarios for dose calculation and optimziation
    pln.multScen = matRad_multScen(ct,'nomScen'); % optimize on the nominal scenario
    
    % Initialisierung der Lungenmodulationsberechnung:
    pln.propHeterogeneity = matRad_HeterogeneityConfig();
    pln.propHeterogeneity.calcHetero = 1; 
    %% Generate Beam Geometry STF
    stf = matRad_generateStf(ct,cst,pln);

    %%
    dij = matRad_calcParticleDose(ct,stf,pln,cst);

    %% Inverse Optimization  for IMPT based on RBE-weighted dose
    resultGUI_homogeneous.patient(i) = matRad_fluenceOptimization(dij,cst,pln);
    %% Redo optimization until V95 = 0.98 für PTV: (hat keine Besserung gebracht!
    % voxIdx = cst{ptvIdx,4}{1,1};
    % roiDose = resultGUI_homogeneous.patient(i).physicalDose(voxIdx);
    % roiDose = sort(roiDose(:), 'descend');
    % V95  = sum(roiDose >= 0.95*cst{ptvIdx,6}{1,1}.parameters{1}/pln.numOfFractions)/numel(roiDose);
    % 
    % while V95<0.98
    %     resultGUI_homogeneous.patient(i) = matRad_fluenceOptimization(dij,cst,pln,resultGUI_homogeneous.patient(i).w);
    % end
    % for n=1:5
    %     resultGUI_homogeneous.patient(i) = matRad_fluenceOptimization(dij,cst,pln,resultGUI_homogeneous.patient(i).w);
    % end
    %% Den Teil auskommentieren für Lungenmodulation weglassen:
    % Anstellen der Heterogenitätsberechnung:
    cst_withLungFlag = pln.propHeterogeneity.cstHeteroAutoassign(cst);
    ct.cube_pmod = pmodCT(i).pmod;
    pln.propHeterogeneity.type = 'local_pmod';
    pln.propHeterogeneity.calcHetero = 1;
    resultGUI_heterogeneous.patient(i) = matRad_calcDoseDirect(ct,stf,pln,cst_withLungFlag,resultGUI_homogeneous.patient(i).w);
    %% Im Folgenden nur noch Darstellung der Ergebnisse und Berechnung von Qualitätsindizes sowie DVH:
    %% Plotten von CT slice axial im Isozentrum:
    isoCenter_mm = matRad_getIsoCenter(cst,ct,0);
    isoCenter_voxel = [round(isoCenter_mm(2)./ct.resolution.y),round(isoCenter_mm(1)./ct.resolution.x),round(isoCenter_mm(3)./ct.resolution.z)];
    x_mm = ((1:size(ct.cubeHU{1,1},2))) * ct.resolution.y;
    y_mm = ((1:size(ct.cubeHU{1,1},1))) * ct.resolution.x;
    filesavefig = "\\david.lse.thm.de\Jessica\MATLAB\Pmod_Patients\CTs\" + patientID +".fig";
    fig = figure;
    imagesc(x_mm, y_mm, ct.cubeHU{1,1}(:,:,isoCenter_voxel(3)))
    title(patientTitle)
    xlabel('x in mm')
    ylabel('y in mm')
    savefig(fig, filesavefig)
    close(fig)
    %% Unterschiede zwischen Plänen: Analyse
    quality_indices = DVH_analysis(resultGUI_homogeneous.patient(i).physicalDose, resultGUI_heterogeneous.patient(i).physicalDose, ct, cst, pln, filedir, "homogenous", "heterogenous", step_size, patientID);
    fieldname = matlab.lang.makeValidName(Patients(i));
    quality_indices_all.(fieldname) = quality_indices;
    % quality_indices.patient(i).name = Patients(i);
    % volume of CTV:
    target = {'CTV'};
    bodyMask = zeros(ct.cubeDim);
    for k = 1:size(cst,1)
        if any(contains(cst{k,2}, target))
            quality_indices_all.(fieldname).CTVvolincm3 = length(cst{k,4}{1,1})*ct.resolution.x*ct.resolution.y*ct.resolution.z*10^(-3); 
            ctvIdx = cst{k,4}{1,1};
        elseif strcmp(cst{k,2}, 'Lungs')
            bodyIdx = cst{k,4}{1,1};
            bodyMask(bodyIdx) = 1;
        elseif any(contains(cst{k,2}, 'PTV'))
            bodyIdx = cst{k,4}{1,1};
            bodyMask(bodyIdx) = 1;
        end
    end
    %% depth of CTV in beam direction:
    step_mm = 1;
    % CTV center
    [xCTV, yCTV, zCTV] = ind2sub(ct.cubeDim, ctvIdx);
    CTVcenter_voxel = [mean(xCTV), mean(yCTV), mean(zCTV)];
    CTVcenter_mm = CTVcenter_voxel .* [ct.resolution.x, ct.resolution.y, ct.resolution.z];
    % helpfunction: check if pos is inside BODY
    isInsideBody = @(pos_mm) ...
        bodyMask( ...
        round(pos_mm(1)/ct.resolution.x), ...
        round(pos_mm(2)/ct.resolution.y), ...
        round(pos_mm(3)/ct.resolution.z) );

    depths_cm = zeros(length(pln.propStf.gantryAngles),1);
    for c = 1:length(pln.propStf.couchAngles)
        for g = 1:length(pln.propStf.gantryAngles)

            % Winkel in Radiant
            phi   = deg2rad(pln.propStf.gantryAngles(g));   % Gantry
            theta = deg2rad(pln.propStf.couchAngles(c));    % Couch

            % Rotation um z (Gantry, clockwise)
            Rz = [ ...
                cos(phi),  -sin(phi), 0; ...
                sin(phi),  cos(phi), 0; ...
                0,         0,        1 ];

            % Rotation um x (Couch, counter-clockwise)
            Ry = [ ...
                1,         0,  0; ...
                0 ,  cos(theta), -sin(theta); ...
                0,  sin(theta), cos(theta) ];

            % Gesamtdrehung
            R = Rz * Ry;

            % Referenzstrahl (AP)
            beam0 = [1; 0;  0];

            % Beam-Richtung im LPS-System
            beamDir = R * beam0;
            beamDir = beamDir / norm(beamDir);
            
            %----Berechnung von CTV Ausdehnung in Strahlrichtung
            CTVbeamdir_voxel = [xCTV, yCTV, zCTV]*beamDir;
            resolution.beamdir = [ct.resolution.x, ct.resolution.y, ct.resolution.z]*beamDir;
            CTVbeamdir_vaus = (max(CTVbeamdir_voxel)-min(CTVbeamdir_voxel));
            CTVbeamdir_cm = (CTVbeamdir_vaus*resolution.beamdir)/10;
            
            %----Berechnung von CTV Ausdehnung senkrecht zur Strahlrichtung
            orth_beamDir = null(beamDir'); % 3x2 Matrix, orthogonal zu d
            senkrecht1 = orth_beamDir(:,1);
            senkrecht2 = orth_beamDir(:,2);
            CTVorthodir1_voxel = [xCTV, yCTV, zCTV]*senkrecht1;
            CTVorthodir2_voxel = [xCTV, yCTV, zCTV]*senkrecht2;
            resolution.orthodir1 = [ct.resolution.x, ct.resolution.y, ct.resolution.z]*senkrecht1;
            resolution.orthodir2 = [ct.resolution.x, ct.resolution.y, ct.resolution.z]*senkrecht2;
            CTVor1 = max(CTVorthodir1_voxel)-min(CTVorthodir1_voxel);
            CTVor2 = max(CTVorthodir2_voxel)-min(CTVorthodir2_voxel);
            CTVor1_cm = CTVor1*resolution.orthodir1/10;
            CTVor2_cm = CTVor2*resolution.orthodir2/10;
            
            % ---- Von CTV nach außen zur Haut gehen ----
            pos_mm = CTVcenter_mm;

            while true
                voxel = round(pos_mm ./ [ct.resolution.x, ct.resolution.y, ct.resolution.z]);

                % Sicherheitscheck (außerhalb CT)
                if any(voxel < 1) || ...
                        voxel(1) > ct.cubeDim(1) || ...
                        voxel(2) > ct.cubeDim(2) || ...
                        voxel(3) > ct.cubeDim(3)
                    break;
                end

                if ~bodyMask(voxel(1), voxel(2), voxel(3))
                    break;
                end

                pos_mm = pos_mm - (beamDir * step_mm)';
            end

            entryPoint_mm = pos_mm;

            % ---- Depth berechnen ----
            depth_mm = norm(CTVcenter_mm - entryPoint_mm);
            depths_cm(c,g) = depth_mm / 10;
            quality_indices_all.(fieldname).gantryAngles(c,g) = pln.propStf.gantryAngles(g);
            quality_indices_all.(fieldname).depths_cm(c,g) = depths_cm(c,g);
            quality_indices_all.(fieldname).CTV_beamdir_cm(c,g) = CTVbeamdir_cm(c,g);
            quality_indices_all.(fieldname).CTV_orthodir1_cm(c,g) = CTVor1_cm(c,g);
            quality_indices_all.(fieldname).CTV_orthodir2_cm(c,g) = CTVor2_cm(c,g);
        end 
        quality_indices_all.(fieldname).couchAngles(c,g) = pln.propStf.couchAngles(c);
    end
    
    
    %% mean Pmod for Beam:
    pmods = ct.cube_pmod(resultGUI_heterogeneous.patient(i).physicalDose~=0);
    quality_indices_all.(fieldname).meanPmodBeam = mean(pmods(pmods~=0)); 
    quality_indices_all.(fieldname).stdPmodBeam = std(pmods(pmods~=0)); 
    quality_indices_all.(fieldname).maxPmodBeam = max(pmods(pmods~=0)); 
    quality_indices_all.(fieldname).minPmodBeam = min(pmods(pmods~=0)); 
    [quality_indices_all.(fieldname).gammaCube,quality_indices_all.(fieldname).gammaPassRate] = compareDose_plots(resultGUI_homogeneous.patient(i).physicalDose, resultGUI_heterogeneous.patient(i).physicalDose,ct,cst,[1 1 0], [], pln, [2,2], 0,'global', filedir_gamma, patientID+"homo_vs_hetero_2%_2mm");

    %% PVH for all patients:
    % Find different positions:
    [idx_lul,~] = ind2sub(size(cst), find(strcmpi(cst,'lung l')));
    [idx_lur,~] = ind2sub(size(cst),find(strcmpi(cst,'lung r')));
    [idx_ctv,~] = ind2sub(size(cst),find(cellfun(@(x) ischar(x) && contains(x,'CTV'), cst)));

    if any(ismember(cell2mat(cst{idx_ctv,4}), cell2mat(cst{idx_lur,4})))
        tumorless_lungwing = "Lung L";
        tumor_lungwing = "Lung R";
    elseif any(ismember(cell2mat(cst{idx_ctv,4}), cell2mat(cst{idx_lul,4})))
        tumorless_lungwing = "Lung R";
        tumor_lungwing = "Lung L";
    else
        matRad_cfg.dispError('No lungwing for Tumor found!');
    end
    %% Define pmod bins
    maxPmod = max(ct.cube_pmod(:));
    PmodBins = 0:step_size_pmod:round(maxPmod,1);

    for j = 1:size(cst,1)
        organName = cst{j,2};
        voxIdx = cst{j,4}{1,1};
        voxIdx_ctv = cst{idx_ctv,4}{1,1};
        voxIdx_without_ctv = voxIdx(~ismember(voxIdx, intersect(voxIdx,voxIdx_ctv)));
        volFrac = (1:numel(voxIdx_without_ctv))/numel(voxIdx_without_ctv);
        % ----- extract pvh values
        if strcmpi(organName,tumor_lungwing)
            roiPmod = ct.cube_pmod(voxIdx_without_ctv);
            roiPmod = sort(roiPmod(:), 'descend');
            if all(roiPmod == 0)
                continue;
            end
            [patient(i).histCounts_pvh, edges_pvh] = histcounts(roiPmod, PmodBins);
            patient(i).histCenters_pvh = (edges_pvh(1:end-1) + edges_pvh(2:end)) / 2;
        end
    end

end

% Cumulative PVH:
set(groot, 'defaultLineLineWidth', 3);         % Standard-Linienstärke
set(groot, 'defaultAxesFontSize', 20);         % Standard-Schriftgröße der Achsen
set(groot, 'defaultAxesLineWidth', 2.5);       % Rahmenbreite der Achsen
set(groot, 'defaultTextFontSize', 20);         % Schriftgröße für Titel und Beschriftungen
set(groot, 'defaultFigureColor', 'w');         % Hintergrundfarbe der Figur (weiß)
set(groot, 'defaultAxesGridLineStyle', '--');  % Stil der Gitternetzlinien
set(groot, 'DefaultTextInterpreter', 'latex');
set(groot, 'DefaultAxesTickLabelInterpreter', 'latex');
set(groot, 'DefaultLegendInterpreter', 'latex');
set(groot, 'DefaultColorbarTickLabelInterpreter', 'latex');
set(groot, 'defaultColorbarFontSize', 20);
set(groot, 'defaultColorbarFontWeight', 'bold');
set(groot, 'defaultColorbarFontAngle', 'normal');
set(groot, 'defaultColorbarLineWidth', 1);
set(groot, 'defaultAxesTitleFontSizeMultiplier', 1);
%% Loop over all patients:
fig = figure;
for i = 1:length(Patients)
    path = fullfile(folder, Patients(i));
    load(char(path))
    token = regexp(Patients(i), '^(Patient_\d+)', 'tokens');
    patientID = token{1}{1};
    patientTitle = strrep(patientID, '_', ' ');
    
    plot(patient(i).histCenters_pvh,patient(i).histCounts_pvh, 'LineStyle','--','DisplayName',string(patientTitle));
    hold on
end
xlabel('$P_{mod}$ in $\mu$m')
ylabel('Frequency')
legend()
% Save figure:
filename = "PmodHist_all_patients_original.fig"; % relative dose difference
fullpath = fullfile(filedir_pvh, filename);
savefig(fig, fullpath)
%
% filename = "Cumulative_DVH_"+ string_comparison +".png";
% fullpath = fullfile(filedir, filename);
% exportgraphics(fig, fullpath, 'Resolution', 1200)
% close(fig);

