

%% run load script for a simple patient 
% load_phantom_run_precalc_modulate_dens
% load 60depth_data_with_modulation_prepared_HLUT_B40s.mat
 load pMod_tables_250_450_800.mat % Resolution 1,5 x 1,5 x 1,5 mm^3
% pmod_dens_800 = readmatrix(['Density_probability_Kilian_cumulative_prob.txt']); %für Resolution: 3,3,3 mm^3 

pmod_dens_800 = interp_poission_pMod800 % acquire which density table should be used
% Stolzenberg implementation (3x3x3 mm^3 resolution), set density zero to very small value:
% pmod_dens_800(1,1) = 0.0000001;

% define the number of randomizations
num_repetitions = 100;
% applying an addition filter for the HU values we are going to modulate
HU_schwelle = [-900 -100];
%modulation.cube{1}(modulation.cubeHU{1}<HU_schwelle(1) | modulation.cubeHU{1}>HU_schwelle(2)) = 0;

% checking for HU values outside the range of [-1000,2995] for the density
% calculation
modulation.cubeHU{2} = ct.cubeHU{1};
modulation.cubeHU{2}(modulation.cubeHU{2} < -1000) = -1000;
modulation.cubeHU{2}(modulation.cubeHU{2} > 2995) = 2995;


%% creating random numbers, masking, reshape 

for i = 1: num_repetitions
    disp(i)
    % Stolzenberg implementation of Kilians model:
    % 1) get a random number from 0 to 1 (out of normal distribution) for
    % all voxels in CT
    nrand = rand(size(ct.cube{1}));
    disp(mean(mean(nrand)))
    modulation.modCube{i}=0.*ct.cube{1};
    % 2) get index of CT according to the corresponding cumulative
    % probability distribution
    for j = 1: size(pmod_dens_800,1)
        % disp(j)
        if j==1
            %case 3x3x3 mm^3:
            %modulation.modCube{i} = modulation.modCube{i} +
            %pmod_dens_800(1,1)*(nrand>0. & nrand<=pmod_dens_800(j,5));
            %case 1,5x1,5x1,5 mm^3:
            modulation.modCube{i} = modulation.modCube{i} + pmod_dens_800(j,1)*(nrand>0. & nrand<=pmod_dens_800(j,5));
        else
            %modulation.modCube{i} = modulation.modCube{i} + pmod_dens_800(j,1)*(nrand>pmod_dens_800(j-1,5) & nrand<=pmod_dens_800(j,5));
            %case 1,5x1,5x1,5 mm^3:
            modulation.modCube{i} = modulation.modCube{i} + pmod_dens_800(j,1)*(nrand>pmod_dens_800(j-1,5) & nrand<=pmod_dens_800(j,5));
        end
    end

    % add original Density values outside the lung
    % modulation.modCube{i}(modulation.cube{1} == 0) = ct.cube{1}(modulation.cube{1} == 0);
    % compute some metrics to check wether the calculations are meaninful
    modulation.metrics.meandata(:,i) = mean(modulation.modCube{i}(modulation.cube{1} == 1));
end

% modulation.metrics.rand_numbers = random_numbers;
modulation.metrics.origCTdens = mean(ct.cube{1}(modulation.cube{1} == 1));
modulation.metrics.num_repetitions = num_repetitions;
% modulation.metrics.poission = interp_poission_pMod;
% modulation.metrics.Pmod = 250;
modulation.metrics.meandens = mean(modulation.metrics.meandata);
