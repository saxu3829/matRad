

%% run load script for a simple patient 
% load_phantom_run_precalc_modulate_dens
% load 60detpht_data_with_modulation_prepared_HLUT_B40s.mat
load pMod_tables_250_450_800.mat

poisson_dens(2,:) =  1:size(poisson_dens,2);
poisson_dens = fliplr(poisson_dens');
interp_poission_pMod = interp_poission_pMod250;
% replace minimum value from 0 density to 1e-6 to make sure TOPAS
% simulation ist run properly! 
% is taken care of in precompute_Densfnc just make sure again!
interp_poission_pMod(interp_poission_pMod == 0) = 0.000001;

%% we continue with modulation instead of ct. 
% cube{2} = all voxels inside lung with its final SPR assigned all
% other voxels = 0
% cube{1} = logical mask with 1 inside lung

% define the number of randomizations
num_repetitions = 100;
% applying an addition filter for the HU values we are going to modulate
HU_schwelle = [-900 -100];
modulation.cube{1}(modulation.cubeHU{1}<HU_schwelle(1) | modulation.cubeHU{1}>HU_schwelle(2)) = 0;

% checking for HU values outside the range of [-1000,2995] for the density
% calculation
modulation.cubeHU{2} = ct.cubeHU{1};
modulation.cubeHU{2}(modulation.cubeHU{2} < -1000) = -1000;
modulation.cubeHU{2}(modulation.cubeHU{2} > 2995) = 2995;


%% modulating the ct here
%%%%%%%%%%%%%%% 1 %%%%%%%%%%%%%%%%%%
% convert HU-Value to density according to MC-Calculation
% Implementation of the following formula:
% needed contants are loaded from file. 
%iv:Ge/Patient/SchneiderHounsfieldUnitSections = 9 -1000 -98 15 23 101 2001 2995 2996 3090
%uv:Ge/Patient/SchneiderDensityOffset = 8 0.00121 1.018 1.03 1.003 1.017 2.201 4.54 1
%uv:Ge/Patient/SchneiderDensityFactor = 8 0.001029700665188 0.000893 0.0 0.001169 0.000592 0.0005 0.0 0.0
% uv:Ge/Patient/SchneiderDensityFactorOffset = 8 1000. 0. 1000. 0. 0. -2000. 0. 0.0
% Formula: Density = (Offset + (Factor*(FactorOffset + HU[-1000,2995] ))) * DensityCorrection
% lookup tables for the density calculation from Schneider et a. 
% load C:\Users\Matth\OneDrive\Matlab\_matRad\lungmodulation\density_schneider.mat
% density correction from Schneider et al., values taken from the TOPAS
% calculation input.
% load C:\Users\Matth\OneDrive\Matlab\_matRad\lungmodulation\DensityCorrection.mat 

%% for testing purpose
% testcube = randi([-900 -100], 5,5,3);
% columnwise_CTcubeHU(:,1)  = testcube(:);

% we use the ct / modulation columnwise to make it faster:
columnwise_CTcubeHU(:,1) = modulation.cubeHU{2}(:);

% use density calculated via HLUT instead of Schneider implementation
columnwise_CTcubeHU(:,2) = ct.cube{1}(:);

% calculation of Density = (Offset + (Factor*(FactorOffset + HU[-1000,2995] ))) * DensityCorrection
% columnwise_CTcubeHU(:,2) = (((interp1(density_schneider(:,1), density_schneider(:,4), columnwise_CTcubeHU(:,1),'previous')...
%     + columnwise_CTcubeHU(:,1) ).*interp1(density_schneider(:,1), density_schneider(:,3), columnwise_CTcubeHU(:,1),'previous') )...
%     + interp1(density_schneider(:,1), density_schneider(:,2), columnwise_CTcubeHU(:,1), 'previous'))...
%     .*DensityCorrection(columnwise_CTcubeHU(:,1)+1001);
% single calcluation steps for check
% Offset(:,1) = interp1(density_schneider(:,1), density_schneider(:,2), HU, 'previous')
% Offset = interp1(density_schneider(:,1), density_schneider(:,2), HU, 'previous')
% Factor(:,1) = interp1(density_schneider(:,1), density_schneider(:,3), HU,'previous')
% Factor = interp1(density_schneider(:,1), density_schneider(:,3), HU,'previous')
% FactorOffset(:,1) = interp1(density_schneider(:,1), density_schneider(:,4), HU,'previous')
% FactorOffset = interp1(density_schneider(:,1), density_schneider(:,4), HU,'previous')
% % Density = (((FactorOffset + HU).*Factor)+Offset).*DensityCorrection;
% ct.dens{1} = reshape(matRad_interp1(hlut(:,1),hlut(:,2),double(ct.cubeHU{1}(:))),ct.cubeDim);
% Density = (((interp1(density_schneider(:,1), density_schneider(:,4), HU,'previous') + HU) .*interp1(density_schneider(:,1), density_schneider(:,3), HU,'previous') )+interp1(density_schneider(:,1), density_schneider(:,2), HU, 'previous'))*DensityCorrection(HU_scaled);

% get index position for look up in poission distr. 
columnwise_CTcubeHU(:,3) = interp1(poisson_dens(:,2),poisson_dens(:,1),columnwise_CTcubeHU(:,2), 'nearest', 'extrap');

% use index for the column and and random number for the row 
% get a a new density from the resulting poisson dens table 
% columnwise_CTcubeHU(:,4) = interp_poission_pMod(sub2ind(size(interp_poission_pMod),randi(1001,size(columnwise_CTcubeHU,1),1),columnwise_CTcubeHU(:,3)));

%% creating random numbers, masking, reshape 
% % add original Density values outside the lung

for i = 1: num_repetitions
    % disp(i)
    % use index for the column and and random number for the row 
    % get a a new density from the resulting poisson dens table
    modulation.modCube{i,1} = reshape(interp_poission_pMod(sub2ind(size(interp_poission_pMod),randi(1001,size(columnwise_CTcubeHU,1),1),columnwise_CTcubeHU(:,3))),ct.cubeDim);
    % add original Density values outside the lung
    modulation.modCube{i,1}(modulation.cube{1} == 0) = ct.cube{1}(modulation.cube{1} == 0);
    % compute some metrics to check wether the calculations are meaninful
    modulation.metrics.meandata(:,i) = mean(modulation.modCube{i,1}(modulation.cube{1} == 1));
end


modulation.metrics.meandens = mean(modulation.metrics.meandata);
% modulation.metrics.rand_numbers = random_numbers;
modulation.metrics.origCTdens = mean(ct.cube{1}(modulation.cube{1} == 1));
modulation.metrics.num_repetitions = num_repetitions;
modulation.metrics.poission = interp_poission_pMod;
modulation.metrics.Pmod = 250;

