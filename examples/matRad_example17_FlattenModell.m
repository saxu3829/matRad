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
% [pmodCT] = pmodBloxx(ct,cst, cluster_size, bin_division, order, edgemethod);
% ct.cube_pmod = pmodCT.cube_pmod;
load('\\david.lse.thm.de\Jessica\MATLAB\Cluster_Algorithm_Marc\Cluster_Robustness.mat')
pmodCT = cs_ct(4);
ct.cube_pmod = pmodCT.pmod;


% Definition of the block size
block_size = [cluster_size, cluster_size, cluster_size];

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
%% 1) Extract left lung out of the original CT

% Create matrix with zeros in the same size as original CT 
pmodCT.cubeHU = zeros(ct.cubeDim);
% Assign HU values of the left lung to the new CT, therefore only the lung
% wing can be seen.
pmodCT.cubeHU(cst{3,4}{1,1}) = ct.cubeHU{1,1}(cst{3,4}{1,1});

% Find maximum expansion of lung voxels
[x, y, z] = ind2sub(size(pmodCT.cubeHU), find(pmodCT.cubeHU~=0)); % Indices of only the left lung in x,y,z
xMin = min(x); xMax = max(x);
yMin = min(y); yMax = max(y);
zMin = min(z); zMax = max(z);

% Crop new CT image to contain only the lung wing
pmodCT.croppedHU = pmodCT.cubeHU(xMin:xMax, yMin:yMax, zMin:zMax);

[x2, y2, z2] = ind2sub(size(pmodCT.croppedHU), find(pmodCT.croppedHU)); % Indices of only the left lung in cropped image

% Show the results of cropping
disp('Originalgröße:');
disp(size(pmodCT.cubeHU));
disp('Zugeschnittene Größe:');
disp(size(pmodCT.croppedHU));
%% Test: plot a slice of the new CT image to verify the correct usage!
image(pmodCT.cubeHU(:,:,50),'CDataMapping','scaled');
xlabel('x in mm');
ylabel('y in mm');
colorbar

%% 2) Cluster voxels of the lung HU values
edgemethod = 'Local'; % Modes for handling the edges; either 'None' or 'Simple' or 'Local' --> LocalX, LocalY, LocalZ possible


% Dimensions of the pmodCT data set (cropped)
[dim1, dim2, dim3] = size(pmodCT.croppedHU);
if strcmp(edgemethod, 'None') | strcmp(edgemethod, 'Simple') 
    % Padding if the original data cannot be divided equally by the block size
    pad1 = mod(block_size(1) - mod(dim1, block_size(1)), block_size(1));
    pad2 = mod(block_size(2) - mod(dim2, block_size(2)), block_size(2));
    pad3 = mod(block_size(3) - mod(dim3, block_size(3)), block_size(3));
    pmodCT.croppedHU_padded = padarray(pmodCT.croppedHU, [pad1, pad2, pad3], 'post');

    % new dimensions
    [new_dim1, new_dim2, new_dim3] = size(pmodCT.croppedHU_padded);

    % reshape the data to the size of blocks
    reshaped_data = reshape(pmodCT.croppedHU_padded, block_size(1), new_dim1 / block_size(1), ...
        block_size(2), new_dim2 / block_size(2), ...
        block_size(3), new_dim3 / block_size(3));

    % output values of each block (pmod values)
    block_pmod = zeros(block_size(1),new_dim1 / block_size(1),block_size(2), new_dim2 / block_size(2), block_size(3),new_dim3 / block_size(3));

    for x = 1:size(reshaped_data, 2)
        for y = 1:size(reshaped_data, 4)
            for z = 1:size(reshaped_data, 6)
                block = reshaped_data(:, x, :, y, :, z); % current block
                if any(block(:) ~= 0) % only process blocks with non-zero elements
                    % disp(block(:))
                    % Find all zero values coming from the edges and handle the
                    % edges in the presented way:
                    % Simple: Look at all neighbouring blocks and find non-zero
                    % values fitting into the missing positions
                    % None: Ignore edges
                    if strcmp(edgemethod, 'Simple')
                        n=1; % Running variable, shows how many block next
                        % to the block1 have been used to fill block1
                        while any(block(:) == 0) % only change block if any value in this block are non-zero
                            k=0; % Second running variable, shows if there is any change in each step
                            % Now come a lot of cases:
                            % Look at all neighboring boxes and it is important
                            % to mind the edges (therefore the ifs)
                            % Loop over all planes (x-y-z- indices)
                            for iz = -n:n
                                for iy = -n:n
                                    for ix = -n:n
                                        % For the planes with already
                                        % investigated values:
                                        if iz <n && iz>-n
                                            if iz<0
                                                if abs(ix) == n
                                                    if ix<0
                                                        if iy<0
                                                            if x+ix> 0  && y+iy> 0 && z+iz> 0
                                                                next_block = reshaped_data(:, x+ix, :, y+iy, :, z+iz);
                                                                k=k+1;
                                                                if any(next_block(block(:) == 0) ~= 0)
                                                                    block(block(:) == 0) = next_block(block(:) == 0);
                                                                else
                                                                    disp('none')
                                                                end
                                                            end
                                                        else
                                                            if x+ix> 0 && y+iy < size(reshaped_data, 4) && z+iz> 0
                                                                next_block = reshaped_data(:, x+ix, :, y+iy, :, z+iz);
                                                                k=k+1;
                                                                if any(next_block(block(:) == 0) ~= 0)
                                                                    block(block(:) == 0) = next_block(block(:) == 0);
                                                                else
                                                                    disp('none')
                                                                end
                                                            end
                                                        end
                                                    elseif ix>0
                                                        if iy<0
                                                            if x+ix < size(reshaped_data, 2) && y+iy> 0 && z+iz > 0
                                                                next_block = reshaped_data(:, x+ix, :, y+iy, :, z+iz);
                                                                k=k+1;
                                                                if any(next_block(block(:) == 0) ~= 0)
                                                                    block(block(:) == 0) = next_block(block(:) == 0);
                                                                else
                                                                    disp('none')
                                                                end
                                                            end
                                                        else
                                                            if x+ix < size(reshaped_data, 2) && y+iy < size(reshaped_data, 4) && z+iz > 0
                                                                next_block = reshaped_data(:, x+ix, :, y+iy, :, z+iz);
                                                                k=k+1;
                                                                if any(next_block(block(:) == 0) ~= 0)
                                                                    block(block(:) == 0) = next_block(block(:) == 0);
                                                                else
                                                                    disp('none')
                                                                end
                                                            end
                                                        end
                                                    end
                                                else
                                                    if abs(iy) == n
                                                        if ix<0
                                                            if iy<0
                                                                if x+ix> 0 && y+iy > 0 && z+iz > 0
                                                                    next_block = reshaped_data(:, x+ix, :, y+iy, :, z+iz);
                                                                    k=k+1;
                                                                    if any(next_block(block(:) == 0) ~= 0)
                                                                        block(block(:) == 0) = next_block(block(:) == 0);
                                                                    else
                                                                        disp('none')
                                                                    end
                                                                end
                                                            elseif iy>0
                                                                if x+ix > 0 && y+iy < size(reshaped_data, 4) && z+iz > 0
                                                                    next_block = reshaped_data(:, x+ix, :, y+iy, :, z+iz);
                                                                    k=k+1;
                                                                    if any(next_block(block(:) == 0) ~= 0)
                                                                        block(block(:) == 0) = next_block(block(:) == 0);
                                                                    else
                                                                        disp('none')
                                                                    end
                                                                end
                                                            end
                                                        else
                                                            if iy<0
                                                                if x+ix < size(reshaped_data, 2) && y+iy > 0 && z+iz > 0
                                                                    next_block = reshaped_data(:, x+ix, :, y+iy, :, z+iz);
                                                                    k=k+1;
                                                                    if any(next_block(block(:) == 0) ~= 0)
                                                                        block(block(:) == 0) = next_block(block(:) == 0);
                                                                    else
                                                                        disp('none')
                                                                    end
                                                                end
                                                            elseif iy>0
                                                                if x+ix < size(reshaped_data, 2) && y+iy < size(reshaped_data, 4) && z+iz > 0
                                                                    next_block = reshaped_data(:, x+ix, :, y+iy, :, z+iz);
                                                                    k=k+1;
                                                                    if any(next_block(block(:) == 0) ~= 0)
                                                                        block(block(:) == 0) = next_block(block(:) == 0);
                                                                    else
                                                                        disp('none')
                                                                    end
                                                                end
                                                            end
                                                        end
                                                    end
                                                end
                                            elseif iz>=0
                                                if abs(ix) == n
                                                    if ix<0
                                                        if iy<0
                                                            if x+ix > 0 && y+iy > 0 && z+iz < size(reshaped_data, 6)
                                                                next_block = reshaped_data(:, x+ix, :, y+iy, :, z+iz);
                                                                k=k+1;
                                                                if any(next_block(block(:) == 0) ~= 0)
                                                                    block(block(:) == 0) = next_block(block(:) == 0);
                                                                else
                                                                    disp('none')
                                                                end
                                                            end
                                                        else
                                                            if x+ix > 0 && y+iy < size(reshaped_data, 4) && z+iz < size(reshaped_data, 6)
                                                                next_block = reshaped_data(:, x+ix, :, y+iy, :, z+iz);
                                                                k=k+1;
                                                                if any(next_block(block(:) == 0) ~= 0)
                                                                    block(block(:) == 0) = next_block(block(:) == 0);
                                                                else
                                                                    disp('none')
                                                                end
                                                            end
                                                        end
                                                    elseif ix>0
                                                        if iy<0
                                                            if x+ix < size(reshaped_data, 2) && y+iy > 0 && z+iz < size(reshaped_data, 6)
                                                                next_block = reshaped_data(:, x+ix, :, y+iy, :, z+iz);
                                                                k=k+1;
                                                                if any(next_block(block(:) == 0) ~= 0)
                                                                    block(block(:) == 0) = next_block(block(:) == 0);
                                                                else
                                                                    disp('none')
                                                                end
                                                            end
                                                        else
                                                            if x+ix < size(reshaped_data, 2) && y+iy < size(reshaped_data, 4) && z+iz < size(reshaped_data, 6)
                                                                next_block = reshaped_data(:, x+ix, :, y+iy, :, z+iz);
                                                                k=k+1;
                                                                if any(next_block(block(:) == 0) ~= 0)
                                                                    block(block(:) == 0) = next_block(block(:) == 0);
                                                                else
                                                                    disp('none')
                                                                end
                                                            end
                                                        end
                                                    end
                                                else
                                                    if abs(iy) == n
                                                        if ix<0
                                                            if iy<0
                                                                if x+ix>0 && y+iy>0 && z+iz < size(reshaped_data, 6)
                                                                    next_block = reshaped_data(:, x+ix, :, y+iy, :, z+iz);
                                                                    k=k+1;
                                                                    if any(next_block(block(:) == 0) ~= 0)
                                                                        block(block(:) == 0) = next_block(block(:) == 0);
                                                                    else
                                                                        disp('none')
                                                                    end
                                                                end
                                                            elseif iy>0
                                                                if x+ix>0 && y+iy < size(reshaped_data, 4) && z+iz < size(reshaped_data, 6)
                                                                    next_block = reshaped_data(:, x+ix, :, y+iy, :, z+iz);
                                                                    k=k+1;
                                                                    if any(next_block(block(:) == 0) ~= 0)
                                                                        block(block(:) == 0) = next_block(block(:) == 0);
                                                                    else
                                                                        disp('none')
                                                                    end
                                                                end
                                                            end
                                                        else
                                                            if iy<0
                                                                if x+ix < size(reshaped_data, 2) && y+iy>0 && z+iz < size(reshaped_data, 6)
                                                                    next_block = reshaped_data(:, x+ix, :, y+iy, :, z+iz);
                                                                    k=k+1;
                                                                    if any(next_block(block(:) == 0) ~= 0)
                                                                        block(block(:) == 0) = next_block(block(:) == 0);
                                                                    else
                                                                        disp('none')
                                                                    end
                                                                end
                                                            elseif iy> 0
                                                                if x+ix < size(reshaped_data, 2) && y+iy < size(reshaped_data, 4) && z+iz < size(reshaped_data, 6)
                                                                    next_block = reshaped_data(:, x+ix, :, y+iy, :, z+iz);
                                                                    k=k+1;
                                                                    if any(next_block(block(:) == 0) ~= 0)
                                                                        block(block(:) == 0) = next_block(block(:) == 0);
                                                                    else
                                                                        disp('none')
                                                                    end
                                                                end
                                                            end
                                                        end
                                                    end
                                                end
                                            end
                                            % New plane z = n
                                        elseif iz == n
                                            if ix<0
                                                if iy<0
                                                    if x+ix > 0 && y+iy > 0 && z+iz < size(reshaped_data, 6)
                                                        next_block = reshaped_data(:, x+ix, :, y+iy, :, z+iz);
                                                        k=k+1;
                                                        if any(next_block(block(:) == 0) ~= 0)
                                                            block(block(:) == 0) = next_block(block(:) == 0);
                                                        else
                                                            disp('none')
                                                        end
                                                    end
                                                else
                                                    if x+ix > 0 && y+iy < size(reshaped_data, 4) && z+iz < size(reshaped_data, 6)
                                                        next_block = reshaped_data(:, x+ix, :, y+iy, :, z+iz);
                                                        k=k+1;
                                                        if any(next_block(block(:) == 0) ~= 0)
                                                            block(block(:) == 0) = next_block(block(:) == 0);
                                                        else
                                                            disp('none')
                                                        end
                                                    end
                                                end
                                            else
                                                if iy<0
                                                    if x+ix < size(reshaped_data, 2) && y+iy > 0 && z+iz < size(reshaped_data, 6)
                                                        next_block = reshaped_data(:, x+ix, :, y+iy, :, z+iz);
                                                        k=k+1;
                                                        if any(next_block(block(:) == 0) ~= 0)
                                                            block(block(:) == 0) = next_block(block(:) == 0);
                                                        else
                                                            disp('none')
                                                        end
                                                    end
                                                else
                                                    if x+ix < size(reshaped_data, 2) && y+iy < size(reshaped_data, 4) && z+iz < size(reshaped_data, 6)
                                                        next_block = reshaped_data(:, x+ix, :, y+iy, :, z+iz);
                                                        k=k+1;
                                                        if any(next_block(block(:) == 0) ~= 0)
                                                            block(block(:) == 0) = next_block(block(:) == 0);
                                                        else
                                                            disp('none')
                                                        end
                                                    end
                                                end
                                            end
                                            % New plane z = -n
                                        elseif iz == -n
                                            if ix<0
                                                if iy<0
                                                    if x+ix> 0  && y+iy> 0 && z+iz> 0
                                                        next_block = reshaped_data(:, x+ix, :, y+iy, :, z+iz);
                                                        k=k+1;
                                                        if any(next_block(block(:) == 0) ~= 0)
                                                            block(block(:) == 0) = next_block(block(:) == 0);
                                                        else
                                                            disp('none')
                                                        end
                                                    end
                                                else
                                                    if x+ix > 0 && y+iy < size(reshaped_data, 4) && z+iz > 0
                                                        next_block = reshaped_data(:, x+ix, :, y+iy, :, z+iz);
                                                        k=k+1;
                                                        if any(next_block(block(:) == 0) ~= 0)
                                                            block(block(:) == 0) = next_block(block(:) == 0);
                                                        else
                                                            disp('none')
                                                        end
                                                    end
                                                end
                                            else
                                                if iy<0
                                                    if x+ix < size(reshaped_data, 2) && y+iy > 0 && z+iz > 0
                                                        next_block = reshaped_data(:, x+ix, :, y+iy, :, z+iz);
                                                        k=k+1;
                                                        if any(next_block(block(:) == 0) ~= 0)
                                                            block(block(:) == 0) = next_block(block(:) == 0);
                                                        else
                                                            disp('none')
                                                        end
                                                    end
                                                else
                                                    if x+ix < size(reshaped_data, 2) && y+iy < size(reshaped_data, 4) && z+iz > 0
                                                        next_block = reshaped_data(:, x+ix, :, y+iy, :, z+iz);
                                                        k=k+1;
                                                        if any(next_block(block(:) == 0) ~= 0)
                                                            block(block(:) == 0) = next_block(block(:) == 0);
                                                        else
                                                            disp('none')
                                                        end
                                                    end
                                                end
                                            end
                                        end
                                    end
                                end
                            end
                            % When changes have been performed in cluster,
                            % coint n to the next level
                            if k>0
                                n=n+1;
                            end
                            % If method does not work wholly on said block:
                            % Pick random numbers out of current voxel to fill
                            % the rest of the voxel
                            % % if k==0
                            % %     values = find(block(:) ~= 0);
                            % %     indices = randperm(numel(values), length(find((block(:) == 0))));
                            % %     block(block(:) == 0) = block(values(indices));
                            % % end
                            % if n>=2 && numel(find(block(:)~=0)) > length(find((block(:) == 0)))
                            %     values = find(block(:) ~= 0);
                            %     indices = randperm(numel(values), length(find((block(:) == 0))));
                            %     block(block(:) == 0) = block(values(indices));
                            % end
                        end
                        disp(num2str(std(block(:))))
                        if any(block(:) == 0)
                            disp('Upsi etwas ist wohl falsch')
                        end
                        %% 3) Perform a gaussian fit
                        % Histogramm values of block for further analysis
                        [counts,edges] = histcounts(block, round(numel(block)/10));
                        centers = (edges(1:end-1) + edges(2:end)) / 2; % Bin center
                        % Change block values to column vector
                        y_f = counts/max(counts); % Histogramm values normalized to maximum counts
                        x_f = centers;
                        % Fit manually using lsqcurvefit
                        gaussEqn = @(p, x) exp(-((x - p(1)).^2) / (2 * p(2)^2));
                        p0 = [max(x_f(y_f==max(y_f))), 20];
                        [p_fit, RESNORM,RESIDUAL,EXITFLAG] = lsqcurvefit(gaussEqn, p0, x_f, y_f);
                        [max1, ind] = max(y_f);
                        y_fit_manual = gaussEqn(p_fit, x_f);
                        % Evaluate goodness-of-fit (e.g., RMSE)
                        rmse_manual = sqrt(mean((y_f - y_fit_manual).^2));
                        pmodCT.rmse(x,y,z) = rmse_manual;
                        disp(['RMSE (Manual Fitting): ', num2str(rmse_manual)]);

                        pmod = real(1.5*(p_fit(2)^2/(-1000*p_fit(1)-p_fit(1)^2))^(1/3)); % pmod in mm
                        block_pmod(:,x, :,y, :,z) = pmod*10^3; % Assign pmod to blocks (in mum)
                        if pmod > 1.0
                            bar(centers,y_f);
                            hold on
                            plot(centers, y_fit_manual, '-.');
                            xlabel('HU values');
                            ylabel('n/$n_{max}$',  'Interpreter', 'latex');
                            title('Cluster size: '+string(cluster_size));
                            legend('data', 'fit')
                            hold off
                            disp(pmod)
                            a=1;
                        end
                        if pmod == 0
                            % Plot histogramm of block
                            bar(centers,y_f);
                            hold on
                            plot(centers, y_fit_manual, '-.');
                            xlabel('HU values');
                            ylabel('n/$n_{max}$',  'Interpreter', 'latex');
                            title('Cluster size: '+string(cluster_size));
                            legend('data', 'fit')
                            hold off
                        end
                        % % Plot histogramm of block
                        % bar(centers,y_f);
                        % hold on
                        % plot(centers, y_fit_manual, '-.');
                        % xlabel('HU values');
                        % ylabel('n/$n_{max}$',  'Interpreter', 'latex');
                        % title('Cluster size: '+string(cluster_size));
                        % legend('data', 'fit')
                        % hold off

                    elseif strcmp(edgemethod, 'None')
                        %  exclude any blocks with zero in them
                        if any(block(:) == 0)
                        else
                            %% 3) Perform a gaussian fit
                            % Histogramm values of block for further analysis
                            [counts,edges] = histcounts(block, 50);
                            centers = (edges(1:end-1) + edges(2:end)) / 2; % Bin center
                            % Change block values to column vector
                            y_f = counts/max(counts); % Histogramm values normalized to maximum counts
                            x_f = centers;

                            % Fit manually using lsqcurvefit
                            gaussEqn = @(p, x) exp(-((x - p(1)).^2) / (2 * p(2)^2));
                            p0 = [x_f(y_f==max(y_f)), 50];
                            [p_fit, RESNORM,RESIDUAL,EXITFLAG] = lsqcurvefit(gaussEqn, p0, x_f, y_f);
                            y_fit_manual = gaussEqn(p_fit, x_f);
                            % Evaluate goodness-of-fit (e.g., RMSE)
                            rmse_manual = sqrt(mean((y_f - y_fit_manual).^2));
                            pmodCT.rmse(x,y,z) = rmse_manual;
                            disp(['RMSE (Manual Fitting): ', num2str(rmse_manual)]);

                            pmod = real(1.5*(p_fit(2)^2/(-1000*p_fit(1)-p_fit(1)^2))^(1/3)); % pmod in mm
                            block_pmod(:,x, :,y, :,z) = pmod*10^3; % Assign pmod to blocks (in mum)
                            disp(pmod)
                            if pmod == 0
                                % Plot histogramm of block
                                bar(centers,y_f);
                                hold on
                                plot(centers, y_fit_manual, '-.');
                                xlabel('HU values');
                                ylabel('n/$n_{max}$',  'Interpreter', 'latex');
                                title('Cluster size: '+string(cluster_size));
                                legend('data', 'fit')
                                hold off
                            end
                            % % Plot histogramm of block
                            % bar(centers,y_f);
                            % hold on
                            % plot(centers, y_fit_manual, '-.');
                            % xlabel('HU values');
                            % ylabel('n/$n_{max}$',  'Interpreter', 'latex');
                            % title('Cluster size: '+string(cluster_size));
                            % legend('data', 'fit')
                            % hold off
                        end
                    end
                end
            end
        end
    end
    % Join all blocks to CT size:
    pmodCT.cropped_pmod_padded = reshape(block_pmod(:), new_dim1, new_dim2, new_dim3);
    % Remove padding from matrix
    pmodCT.cropped_pmod = pmodCT.cropped_pmod_padded(1:end-pad1, 1:end-pad2, 1:end-pad3);

%% Local calculation of pmod for each point in the lung
% Steps: 
% 1) Loop over all Voxels in CT and look only at voxels where the value %%
%% anpassen!!
% is non-zero
% 1) Loop over all voxels in CT that are non-zero
% 2) For each value seach the neighbouring voxels for non-zero values and assign them to a matrix
% with the cluster size until the matrix is filled
% 3) Perform cluster analysis for this matrix
% 4) Continue onto next value and repeat 1) to 3) until the end of CT
elseif strcmp(edgemethod, 'Local')
    for x = 1:size(pmodCT.croppedHU, 1)
        disp('x: '+string(x))
        for y = 1:size(pmodCT.croppedHU, 2)
            disp('y: '+string(y))
            for z = 1:size(pmodCT.croppedHU, 3)
                disp('z: '+string(z))
                tic
                if pmodCT.croppedHU(x,y,z) ~=0
                    HU_matrix = zeros(cluster_size, cluster_size, cluster_size);
                    HU_matrix(1) = pmodCT.croppedHU(x,y,z);
                    index = 2;
                    maximum_index = cluster_size*cluster_size*cluster_size;
                    n=1; % Running variable, shows how many neighboring values 
                    % have been addressed already 
                    while any(HU_matrix(:) == 0)
                        % Now come a lot of cases:
                        % Look at all neighboring boxes and it is important
                        % to mind the edges (therefore the ifs)
                        k=0; % Second running variable, shows if there is any change in each step
                        % Now come a lot of cases:
                        % Look at all neighboring boxes and it is important
                        % to mind the edges (therefore the ifs)
                        % Loop over all planes (x-y-z- indices)
                        for iz = -n:n
                            for iy = -n:n
                                for ix = -n:n
                                    % For the planes with already
                                    % investigated values:
                                    if iz <n && iz>-n
                                        if iz<0
                                            if abs(ix) == n
                                                if ix<0
                                                    if iy<0
                                                        if x+ix> 0  && y+iy> 0 && z+iz> 0
                                                            value = pmodCT.croppedHU(x+ix, y+iy, z+iz);
                                                            k=k+1;
                                                            if value ~= 0 && index<=maximum_index
                                                                HU_matrix(index) = value;
                                                                index = index + 1;
                                                            end
                                                        end
                                                    else
                                                        if x+ix> 0 && y+iy < size(pmodCT.croppedHU, 2) && z+iz> 0
                                                            value = pmodCT.croppedHU(x+ix, y+iy, z+iz);
                                                            k=k+1;
                                                            if value ~= 0 && index<=maximum_index
                                                                HU_matrix(index) = value;
                                                                index = index + 1;
                                                            end
                                                        end
                                                    end
                                                elseif ix>0
                                                    if iy<0
                                                        if x+ix < size(pmodCT.croppedHU, 1) && y+iy> 0 && z+iz > 0
                                                            value = pmodCT.croppedHU(x+ix, y+iy, z+iz);
                                                            k=k+1;
                                                            if value ~= 0 && index<=maximum_index
                                                                HU_matrix(index) = value;
                                                                index = index + 1;
                                                            end
                                                        end
                                                    else
                                                        if x+ix < size(pmodCT.croppedHU, 1) && y+iy < size(pmodCT.croppedHU, 2) && z+iz > 0
                                                            value = pmodCT.croppedHU(x+ix, y+iy, z+iz);
                                                            k=k+1;
                                                            if value ~= 0 && index<=maximum_index
                                                                HU_matrix(index) = value;
                                                                index = index + 1;
                                                            end
                                                        end
                                                    end
                                                end
                                            else
                                                if abs(iy) == n
                                                    if ix<0
                                                        if iy<0
                                                            if x+ix> 0 && y+iy > 0 && z+iz > 0
                                                                value = pmodCT.croppedHU(x+ix, y+iy, z+iz);
                                                                k=k+1;
                                                                if value ~= 0 && index<=maximum_index
                                                                    HU_matrix(index) = value;
                                                                    index = index + 1;
                                                                end
                                                            end
                                                        elseif iy>0
                                                            if x+ix > 0 && y+iy < size(pmodCT.croppedHU, 2) && z+iz > 0
                                                                value = pmodCT.croppedHU(x+ix, y+iy, z+iz);
                                                                k=k+1;
                                                                if value ~= 0 && index<=maximum_index
                                                                    HU_matrix(index) = value;
                                                                    index = index + 1;
                                                                end
                                                            end
                                                        end
                                                    else
                                                        if iy<0
                                                            if x+ix < size(pmodCT.croppedHU, 1) && y+iy > 0 && z+iz > 0
                                                                value = pmodCT.croppedHU(x+ix, y+iy, z+iz);
                                                                k=k+1;
                                                                if value ~= 0 && index<=maximum_index
                                                                    HU_matrix(index) = value;
                                                                    index = index + 1;
                                                                end
                                                            end
                                                        elseif iy>0
                                                            if x+ix < size(pmodCT.croppedHU, 1) && y+iy < size(pmodCT.croppedHU, 2) && z+iz > 0
                                                                value = pmodCT.croppedHU(x+ix, y+iy, z+iz);
                                                                k=k+1;
                                                                if value ~= 0 && index<=maximum_index
                                                                    HU_matrix(index) = value;
                                                                    index = index + 1;
                                                                end
                                                            end
                                                        end
                                                    end
                                                end
                                            end
                                        elseif iz>=0
                                            if abs(ix) == n
                                                if ix<0
                                                    if iy<0
                                                        if x+ix > 0 && y+iy > 0 && z+iz < size(pmodCT.croppedHU, 3)
                                                            value = pmodCT.croppedHU(x+ix, y+iy, z+iz);
                                                            k=k+1;
                                                            if value ~= 0 && index<=maximum_index
                                                                HU_matrix(index) = value;
                                                                index = index + 1;
                                                            end
                                                        end
                                                    else
                                                        if x+ix > 0 && y+iy < size(pmodCT.croppedHU, 2) && z+iz < size(pmodCT.croppedHU, 3)
                                                            value = pmodCT.croppedHU(x+ix, y+iy, z+iz);
                                                            k=k+1;
                                                            if value ~= 0 && index<=maximum_index
                                                                HU_matrix(index) = value;
                                                                index = index + 1;
                                                            end
                                                        end
                                                    end
                                                elseif ix>0
                                                    if iy<0
                                                        if x+ix < size(pmodCT.croppedHU, 1) && y+iy > 0 && z+iz < size(pmodCT.croppedHU, 3)
                                                            value = pmodCT.croppedHU(x+ix, y+iy, z+iz);
                                                            k=k+1;
                                                            if value ~= 0 && index<=maximum_index
                                                                HU_matrix(index) = value;
                                                                index = index + 1;
                                                            end
                                                        end
                                                    else
                                                        if x+ix < size(pmodCT.croppedHU, 1) && y+iy < size(pmodCT.croppedHU, 2) && z+iz < size(pmodCT.croppedHU, 3)
                                                            value = pmodCT.croppedHU(x+ix, y+iy, z+iz);
                                                            k=k+1;
                                                            if value ~= 0 && index<=maximum_index
                                                                HU_matrix(index) = value;
                                                                index = index + 1;
                                                            end
                                                        end
                                                    end
                                                end
                                            else
                                                if abs(iy) == n
                                                    if ix<0
                                                        if iy<0
                                                            if x+ix>0 && y+iy>0 && z+iz < size(pmodCT.croppedHU, 3)
                                                                value = pmodCT.croppedHU(x+ix, y+iy, z+iz);
                                                                k=k+1;
                                                                if value ~= 0 && index<=maximum_index
                                                                    HU_matrix(index) = value;
                                                                    index = index + 1;
                                                                end
                                                            end
                                                        elseif iy>0
                                                            if x+ix>0 && y+iy < size(pmodCT.croppedHU, 2) && z+iz < size(pmodCT.croppedHU, 3)
                                                                value = pmodCT.croppedHU(x+ix, y+iy, z+iz);
                                                                k=k+1;
                                                                if value ~= 0 && index<=maximum_index
                                                                    HU_matrix(index) = value;
                                                                    index = index + 1;
                                                                end
                                                            end
                                                        end
                                                    else
                                                        if iy<0
                                                            if x+ix < size(pmodCT.croppedHU, 1) && y+iy>0 && z+iz < size(pmodCT.croppedHU, 3)
                                                                value = pmodCT.croppedHU(x+ix, y+iy, z+iz);
                                                                k=k+1;
                                                                if value ~= 0 && index<=maximum_index
                                                                    HU_matrix(index) = value;
                                                                    index = index + 1;
                                                                end
                                                            end
                                                        elseif iy> 0
                                                            if x+ix < size(pmodCT.croppedHU, 1) && y+iy < size(pmodCT.croppedHU, 2) && z+iz < size(pmodCT.croppedHU, 3)
                                                                value = pmodCT.croppedHU(x+ix, y+iy, z+iz);
                                                                k=k+1;
                                                                if value ~= 0 && index<=maximum_index
                                                                    HU_matrix(index) = value;
                                                                    index = index + 1;
                                                                end
                                                            end
                                                        end
                                                    end
                                                end
                                            end
                                        end
                                        % New plane z = n
                                    elseif iz == n
                                        if ix<0
                                            if iy<0
                                                if x+ix > 0 && y+iy > 0 && z+iz < size(pmodCT.croppedHU, 3)
                                                    value = pmodCT.croppedHU(x+ix, y+iy, z+iz);
                                                    k=k+1;
                                                    if value ~= 0 && index<=maximum_index
                                                        HU_matrix(index) = value;
                                                        index = index + 1;
                                                    end
                                                end
                                            else
                                                if x+ix > 0 && y+iy < size(pmodCT.croppedHU, 2) && z+iz < size(pmodCT.croppedHU, 3)
                                                    value = pmodCT.croppedHU(x+ix, y+iy, z+iz);
                                                    k=k+1;
                                                    if value ~= 0 && index<=maximum_index
                                                        HU_matrix(index) = value;
                                                        index = index + 1;
                                                    end
                                                end
                                            end
                                        else
                                            if iy<0
                                                if x+ix < size(pmodCT.croppedHU, 1) && y+iy > 0 && z+iz < size(pmodCT.croppedHU, 3)
                                                    value = pmodCT.croppedHU(x+ix, y+iy, z+iz);
                                                    k=k+1;
                                                    if value ~= 0 && index<=maximum_index
                                                        HU_matrix(index) = value;
                                                        index = index + 1;
                                                    end
                                                end
                                            else
                                                if x+ix < size(pmodCT.croppedHU, 1) && y+iy < size(pmodCT.croppedHU, 2) && z+iz < size(pmodCT.croppedHU, 3)
                                                    value = pmodCT.croppedHU(x+ix, y+iy, z+iz);
                                                    k=k+1;
                                                    if value ~= 0 && index<=maximum_index
                                                        HU_matrix(index) = value;
                                                        index = index + 1;
                                                    end
                                                end
                                            end
                                        end
                                        % New plane z = -n
                                    elseif iz == -n
                                        if ix<0
                                            if iy<0
                                                if x+ix> 0  && y+iy> 0 && z+iz> 0
                                                    value = pmodCT.croppedHU(x+ix, y+iy, z+iz);
                                                    k=k+1;
                                                    if value ~= 0 && index<=maximum_index
                                                        HU_matrix(index) = value;
                                                        index = index + 1;
                                                    end
                                                end
                                            else
                                                if x+ix > 0 && y+iy < size(pmodCT.croppedHU, 2) && z+iz > 0
                                                    value = pmodCT.croppedHU(x+ix, y+iy, z+iz);
                                                    k=k+1;
                                                    if value ~= 0 && index<=maximum_index
                                                        HU_matrix(index) = value;
                                                        index = index + 1;
                                                    end
                                                end
                                            end
                                        else
                                            if iy<0
                                                if x+ix < size(pmodCT.croppedHU, 1) && y+iy > 0 && z+iz > 0
                                                    value = pmodCT.croppedHU(x+ix, y+iy, z+iz);
                                                    k=k+1;
                                                    if value ~= 0 && index<=maximum_index
                                                        HU_matrix(index) = value;
                                                        index = index + 1;
                                                    end
                                                end
                                            else
                                                if x+ix < size(pmodCT.croppedHU, 1) && y+iy < size(pmodCT.croppedHU, 2) && z+iz > 0
                                                    value = pmodCT.croppedHU(x+ix, y+iy, z+iz);
                                                    k=k+1;
                                                    if value ~= 0 && index<=maximum_index
                                                        HU_matrix(index) = value;
                                                        index = index + 1;
                                                    end
                                                end
                                            end
                                        end
                                    end
                                end
                            end
                        end
                        % When changes have been performed in cluster,
                        % coint n to the next level
                        if k>0
                            n=n+1;
                        end
                        % If method does not work wholly on said block:
                        % Pick random numbers out of current voxel to fill
                        % the rest of the voxel
                        if k==0 | n>21
                            values = find(HU_matrix(:) ~= 0);
                            indices = randperm(numel(values), length(find((HU_matrix(:) == 0))));
                            HU_matrix(HU_matrix(:) == 0) = HU_matrix(values(indices));
                        end
                        % if n>=2 && numel(find(block(:)~=0)) > length(find((block(:) == 0)))
                        %     values = find(block(:) ~= 0);
                        %     indices = randperm(numel(values), length(find((block(:) == 0))));
                        %     block(block(:) == 0) = block(values(indices));
                        % end
                    end
                    if index<maximum_index
                        disp('FEHLER')
                    end    
                    disp(num2str(std(HU_matrix(:))))
                    if any(HU_matrix(:) == 0)
                        disp('Upsi etwas ist wohl falsch')
                    end
                    %% 3) Perform a gaussian fit
                    % Histogramm values of block for further analysis
                    [counts,edges] = histcounts(HU_matrix, round(numel(HU_matrix)/10));
                    centers = (edges(1:end-1) + edges(2:end)) / 2; % Bin center
                    % Change block values to column vector
                    y_f = counts/max(counts); % Histogramm values normalized to maximum counts
                    x_f = centers;
                    % Fit manually using lsqcurvefit
                    gaussEqn = @(p, x) exp(-((x - p(1)).^2) / (2 * p(2)^2));
                    p0 = [max(x_f(y_f==max(y_f))), 20];
                    [p_fit, RESNORM,RESIDUAL,EXITFLAG] = lsqcurvefit(gaussEqn, p0, x_f, y_f);
                    [max1, ind] = max(y_f);
                    y_fit_manual = gaussEqn(p_fit, x_f);
                    % Evaluate goodness-of-fit (e.g., RMSE)
                    rmse_manual = sqrt(mean((y_f - y_fit_manual).^2));
                    pmodCT.rmse(x,y,z) = rmse_manual;
                    disp(['RMSE (Manual Fitting): ', num2str(rmse_manual)]);

                    pmod = real(1.5*(p_fit(2)^2/(-1000*p_fit(1)-p_fit(1)^2))^(1/3)); % pmod in mm
                    pmodCT.cropped_pmod(x,y,z) = pmod*10^3; % Assign pmod to blocks (in mum)
                    if pmod > 1.0
                        bar(centers,y_f);
                        hold on
                        plot(centers, y_fit_manual, '-.');
                        xlabel('HU values');
                        ylabel('n/$n_{max}$',  'Interpreter', 'latex');
                        title('Cluster size: '+string(cluster_size));
                        legend('data', 'fit')
                        hold off
                        disp(pmod)
                        a=1;
                    end
                    if pmod == 0
                        % Plot histogramm of block
                        bar(centers,y_f);
                        hold on
                        plot(centers, y_fit_manual, '-.');
                        xlabel('HU values');
                        ylabel('n/$n_{max}$',  'Interpreter', 'latex');
                        title('Cluster size: '+string(cluster_size));
                        legend('data', 'fit')
                        hold off
                    end
                    % % Plot histogramm of block
                    % bar(centers,y_f);
                    % hold on
                    % plot(centers, y_fit_manual, '-.');
                    % xlabel('HU values');
                    % ylabel('n/$n_{max}$',  'Interpreter', 'latex');
                    % title('Cluster size: '+string(cluster_size));
                    % legend('data', 'fit')
                    % hold off
                end
                toc
            end
        end
    end
end


%% Now put the pmods in place of the previous HU values: error needs to be solved!
% Find missing dimensions of cropped array(pad zeros to cropped ct):
cube_pmod = pmodCT.cropped_pmod;
if ct.cubeDim(1) > dim1
    cube_pmod = padarray(cube_pmod, [xMin-1 0 0], 'pre');
    cube_pmod = padarray(cube_pmod, [ct.cubeDim(1)-xMax 0 0], 'post');
end
if ct.cubeDim(2) > dim2
    cube_pmod = padarray(cube_pmod, [0 yMin-1 0], 'pre');
    cube_pmod = padarray(cube_pmod, [0 ct.cubeDim(2)-yMax 0], 'post');
end
if ct.cubeDim(3) > dim3
    cube_pmod = padarray(cube_pmod, [0 0 zMin-1], 'pre');
    cube_pmod = padarray(cube_pmod, [0 0 ct.cubeDim(3)-zMax], 'post');
end

% HU range for the analysis of lung structures(all else no pmod determined):
HU_range = [-900 -100];
pmodCT.cubeHU(pmodCT.cubeHU<HU_range(1) | pmodCT.cubeHU>HU_range(2)) = 0;

% Assign pmod values to CT:
pmodCT.cube_pmod = pmodCT.cubeHU;

pmodCT.cube_pmod(pmodCT.cubeHU~=0) = real(cube_pmod(pmodCT.cubeHU~=0));


% Plot resulting Pmod slice:
figure
image(pmodCT.cube_pmod(:,:,50),'CDataMapping','scaled');
xlabel('x in mm');
ylabel('y in mm');
c = colorbar;   % Colorbar hinzufügen
c.Label.String = 'Pmod in \mum'; % Label setzen
set(gca, 'ColorScale', 'log');
% Plot histogram of pmod values in CT that are non-zero:

figure
histogram(pmodCT.cube_pmod(pmodCT.cube_pmod~=0));
set(gca, 'YScale', 'log');
xlabel('Pmod in $\mu m$',  'Interpreter', 'latex');
ylabel('n');
% Plot histogram of RMSE values in CT that are non-zero:
figure
histogram(pmodCT.rmse(pmodCT.rmse~=0),10);
set(gca, 'YScale', 'log');
xlabel('RMSE',  'Interpreter', 'latex');
ylabel('n');

toc