function quality_indices = DVH_analysis(cube1, cube2, ct, cst, pln, filedir, planA, planB, step_size, Patientname)
% Calculation of dose quality metrics and quality indices to compare plans
% ( robustness analysis, inspiration by Kim et al. (2025) table 2 )
%
% call
%   [quality_indices1, quality_indices2] = DVH_analysis(cube1, cube2, ct, cst, pln, criteria, filedir, string_comparison)
%
% input
%   cube1:         dose cube 1 as an M x N x O array
%   cube2:         dose cube 2 as an M x N x O array
%   ct:            ct struct with ct cube
%   cst:           list of interesting volumes inside the patient as matRad
%                  struct (optional, does not calculate gamma and DVH)
%   pln:            to know the number of fractions used to calculate the
%   plan
%   filedir:       directory to save plots
%   planA:         string which specifies cube1 ( default: plan1 )
%   planB:         string which specifies cube2 ( default: plan2 ) 
%   step_size (optionial):     step size of dvh histogramm (dose in Gy) 
%                              default: 0.1 Gy
%
% output
%
%   quality_indices1/2: calculated quality indices (homogeneity_index,
%   conformity_index, CTV D95%, etc. for cube1/2
%
% References:
%   [1]  https://doi.org/10.1016/j.clon.2025.103938
%   [2]  Homogeneity Index: 
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2015 the matRad development team.
%
% This file is part of the matRad project. It is subject to the license
% terms in the LICENSE file found in the top-level directory of this
% distribution and at https://github.com/e0404/matRad/LICENSES.txt. No part
% of the matRad project, including this file, may be copied, modified,
% propagated, or distributed except according to the terms contained in the
% LICENSE file.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

matRad_cfg = MatRad_Config.instance();
if nargin < 10
    step_size = 0.1;
end
if nargin < 7 || isempty(planA)
    planA = "plan1";
end
if nargin < 8 || isempty(planB)
    planB = "plan2";
end

%% Consistency check
if ~isequal(size(cube1), size(cube2))
    matRad_cfg.dispError('Dose cubes must be the same size\n');
end

if ~exist('cst', 'var') || isempty(cst)
    skip = 1;
    skip_pvh = 0;
    matRad_cfg.dispWarning('No CST provided — skipping DVH');
elseif ~isfield(ct, 'cube_pmod') || isempty(ct.cube_pmod)
    matRad_cfg.dispWarning('No PmodMap provided — skipping PVH');
    skip = 0;
    skip_pvh = 1;
else
    skip = 0;
    skip_pvh = 0;
end
if ~exist('pln', 'var') || isempty(pln)
    pln = [];
end

if ~skip
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

    %% --- Identify reference dose from CTV/PTV ---
    target = {'CTV','PTV'};
    for i = 1:size(cst,1)
        if any(contains(cst{i,2}, target)) & ~isempty(cst{i,6})
            ct.planDose = cst{i,6}{1,1}.parameters{1}/pln.numOfFractions;
        end
    end
    lung = {'Lung'};
    marker = {'BODY','Isocenter Marker'};

    %% --- Define dose bins ---
    maxDose1 = max(cube1(:));
    maxDose2 = max(cube2(:));
    doseBins{1} = 0:step_size:round(maxDose1,1);
    doseBins{2} = 0:step_size:round(maxDose2,1);

    planNames = {planA, planB};
    doseCubes = {cube1,cube2};
   
    %% --- Loop over structures ---
    for j = 1:size(cst,1)
        organName = cst{j,2};
        if any(contains(organName, marker))
            if any(contains(organName, 'BODY'))
                organField = matlab.lang.makeValidName(organName);
                voxIdx = cst{j,4}{1,1};
                volFrac = (1:numel(voxIdx))/numel(voxIdx);
                quality_indices.(organField).organ = organName;
                for p = 1:2
                    % --- extract dose values ---
                    roiDose = doseCubes{p}(voxIdx);
                    roiDose = sort(roiDose(:), 'descend');
                    if all(roiDose == 0)
                        continue;
                    end
                    Q.V95  = sum(roiDose >= 0.95*ct.planDose)/numel(roiDose);

                    % --- assign ---
                    quality_indices.(organField).(planNames{p}) = Q;
                end
            else
                continue;
            end
        end
        organField = matlab.lang.makeValidName(organName);
        voxIdx = cst{j,4}{1,1};
        volFrac = (1:numel(voxIdx))/numel(voxIdx);
        quality_indices.(organField).organ = organName;

        for p = 1:2
            % --- extract dose values ---
            roiDose = doseCubes{p}(voxIdx);
            roiDose = sort(roiDose(:), 'descend');
            if all(roiDose == 0)
                continue;
            end
            
            % --- DVH bins (optional, if you want to store) ---
            [histCounts, edges] = histcounts(roiDose, doseBins{p});
            histCenters = (edges(1:end-1) + edges(2:end)) / 2;
            dvh.(organField).(planNames{p}).dose = roiDose;
            dvh.(organField).(planNames{p}).centers = histCenters;
            dvh.(organField).(planNames{p}).counts = histCounts;

            % --- Dose metrics ---
            Q.Dmax = max(roiDose);
            Q.Dmin = min(roiDose);
            Q.Dmean = mean(roiDose);
            Q.D98  = roiDose(find(volFrac >= 0.98,1));
            Q.D95  = roiDose(find(volFrac >= 0.95,1));
            Q.D50  = roiDose(find(volFrac >= 0.5,1));
            Q.D5   = roiDose(find(volFrac >= 0.05,1));
            Q.D2   = roiDose(find(volFrac >= 0.02,1));

            % --- Volume metrics ---
            Q.V95  = sum(roiDose >= 0.95*ct.planDose)/numel(roiDose);
            Q.V100 = sum(roiDose >= ct.planDose)/numel(roiDose);

            % --- Homogeneity & OAR metrics ---
            if any(contains(organName, target))
                Q.HI = (Q.D2 - Q.D98) / Q.D50;
            elseif any(contains(organName, lung)) && ~strcmp(organName, tumorless_lungwing)
                Q.V20Gy = sum(roiDose .* pln.numOfFractions >= 20) / numel(roiDose);
            end

            % --- assign ---
            quality_indices.(organField).(planNames{p}) = Q;
        end

        % --- compute deltas ---
        if isfield(quality_indices.(organField), planA) && isfield(quality_indices.(organField), planB)
            fieldsQ = fieldnames(quality_indices.(organField).(planA));
            for f = 1:numel(fieldsQ)
                fName = fieldsQ{f};
                valA = quality_indices.(organField).(planA).(fName);
                valB = quality_indices.(organField).(planB).(fName);
                if isnumeric(valA) && isnumeric(valB) && valA~=0
                    quality_indices.(organField).Delta.(fName) = (abs(1 - valB / valA))*100; % in %
                elseif isnumeric(valA) && isnumeric(valB)
                    quality_indices.(organField).Delta.(fName) = valB;
                end
            end
        end
    end

    %% --- Conformity Index & Undercoverage ---
    ptv = {'PTV'}; body = {'BODY'};
    try
        ptvIdx = find(contains(cst(:,2), ptv));
        bodyIdx = find(contains(cst(:,2), body));
        ptvVox = numel(cst{ptvIdx,4}{1,1});
        bodyVox = numel(cst{bodyIdx,4}{1,1});
        PTVField = matlab.lang.makeValidName(cst{ptvIdx,2});
        BODYField = matlab.lang.makeValidName(cst{bodyIdx,2});

        for p = 1:2
            planN = planNames{p};
            quality_indices.(PTVField).(planN).CI_PTV = ...
                quality_indices.(BODYField).(planN).V95*bodyVox / ptvVox;
            quality_indices.(PTVField).(planN).Vunder = ...
                (1 - quality_indices.(PTVField).(planN).V95) * 100;
        end

        quality_indices.(PTVField).Delta.CI_PTV = ...
            (abs(1 - quality_indices.(PTVField).(planB).CI_PTV / quality_indices.(PTVField).(planA).CI_PTV))*100; % in %
        quality_indices.(PTVField).Delta.Vunder = ...
            abs(quality_indices.(PTVField).(planB).Vunder - quality_indices.(PTVField).(planA).Vunder); % in %
    catch
        warning('PTV/CTV not found for Conformity Index calculation.');
    end
    
    %% Plot DVHs
    set(groot, 'defaultLineLineWidth', 3);         % Standard-Linienstärke
    set(groot, 'defaultAxesFontSize', 20);         % Standard-Schriftgröße der Achsen
    set(groot, 'defaultAxesLineWidth', 1.5);       % Rahmenbreite der Achsen
    set(groot, 'defaultTextFontSize', 15);         % Schriftgröße für Titel und Beschriftungen
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

    % % Differential DVH:
    % fig = figure('Units','normalized','Position',[0. 0. 1. 1.]);
    % for j = 1:length(cst(:,1))
    %     if all(cellfun(@(teststr) isempty(strfind(cst{j,2},teststr)), marker)) % exclude isocenter marker
    %         if any(dvh.dosevoxel1{j}~=0) || any(dvh.dosevoxel2{j}~=0)
    %             plot(dvh.histcenters1{j},dvh.histcounts1{j}.*100./sum(dvh.histcounts1{j}), 'Color', cst{j,5}.visibleColor,'LineStyle','--','DisplayName',string(cst{j,2})+"_Cube1");
    %             hold on
    %             %plot(dvh.dvhistcounts2{j},dvh.volume{j}, 'Color', cst{j,5}.visibleColor,'LineStyle',':','DisplayName',string(cst{j,2})+"_Cube2");
    %             plot(dvh.histcenters2{j},dvh.histcounts2{j}.*100./sum(dvh.histcounts2{j}), 'Color', cst{j,5}.visibleColor,'LineStyle',':','DisplayName',string(cst{j,2})+"_Cube2");
    %         end
    %     end
    % end
    % xlabel('Dose in Gy')
    % ylabel('Volume in %')
    % legend()
    % % Save figure:
    % filename = "Differential_DVH_"+ string_comparison +".fig"; % relative dose difference
    % fullpath = fullfile(filedir, filename);
    % % savefig(fig, fullpath)
    %
    % % filename = "Cumulative_DVH_"+ string_comparison +".png";
    % % fullpath = fullfile(filedir, filename);
    % % exportgraphics(fig, fullpath, 'Resolution', 1200)
    % % close(fig);                  % Figure schließen

    % Cumulative DVH:


    fig = figure('Units','normalized','Position',[0. 0. 1. 1.]);
    for j = 1:length(cst(:,1))
        organName = cst{j,2};
        organField = matlab.lang.makeValidName(organName);
        if any(contains(organName, marker))
            continue;
        end
        if any(doseCubes{1}(cst{j,4}{1,1})~=0) || any(doseCubes{2}(cst{j,4}{1,1})~=0)
            if isequal(cst{j,5}.visibleColor, [1 1 1]) % schließt weiß als Farbe aus und setzt es zu schwarz
                plot(dvh.(organField).(planNames{1}).centers,cumsum(dvh.(organField).(planNames{1}).counts.*100./sum(dvh.(organField).(planNames{1}).counts), 'reverse'), 'Color', [0 0 0],'LineStyle','--','DisplayName',string(organField)+"_"+string(planA));
                hold on
                plot(dvh.(organField).(planNames{2}).centers,cumsum(dvh.(organField).(planNames{2}).counts.*100./sum(dvh.(organField).(planNames{2}).counts), 'reverse'), 'Color', [0 0 0],'LineStyle',':','DisplayName',string(organField)+"_"+string(planB));
            else
                plot(dvh.(organField).(planNames{1}).centers,cumsum(dvh.(organField).(planNames{1}).counts.*100./sum(dvh.(organField).(planNames{1}).counts), 'reverse'), 'Color', cst{j,5}.visibleColor,'LineStyle','--','DisplayName',string(organField)+"_"+string(planA));
                hold on
                plot(dvh.(organField).(planNames{2}).centers,cumsum(dvh.(organField).(planNames{2}).counts.*100./sum(dvh.(organField).(planNames{2}).counts), 'reverse'), 'Color', cst{j,5}.visibleColor,'LineStyle',':','DisplayName',string(organField)+"_"+string(planB));
            end
        end
    end
    xlabel('Dose in Gy')
    ylabel('Volume in %')
    title(Patientname)
    legend()
    % Save figure:
    filename = "Cumulative_DVH_"+ string(Patientname) + planA + "vs" + planB +".fig"; % relative dose difference
    fullpath = fullfile(filedir, filename);
    savefig(fig, fullpath)
    %
    % filename = "Cumulative_DVH_"+ string_comparison +".png";
    % fullpath = fullfile(filedir, filename);
    % exportgraphics(fig, fullpath, 'Resolution', 1200)
    close(fig);                  % Figure schließen

    matRad_cfg.dispInfo('Done!\n');
end
end        

