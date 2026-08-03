function bixel = matRad_calcParticleDoseBixel(radDepths, radialDist_sq, sigmaIni_sq, baseData, heteroCorrDepths, propHeterogeneity , vTissueIndex, pmod)
% matRad visualization of two-dimensional dose distributions
% on ct including segmentation
%
% call
%   dose = matRad_calcParticleDoseBixel(radDepths, radialDist_sq, sigmaIni_sq, baseData)
%
% input
%   radDepths:      radiological depths
%   radialDist_sq:  squared radial distance in BEV from central ray
%   sigmaIni_sq:    initial Gaussian sigma^2 of beam at patient surface
%   baseData:       base data required for particle dose calculation
%
% output
%   dose:   particle dose at specified locations as linear vector
%
% References
%   [1] http://iopscience.iop.org/0031-9155/41/8/005
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

% Instance of MatRad_Config class
matRad_cfg = MatRad_Config.instance();

% skip heterogeneity correction for other functions
if nargin < 5
    heteroCorrDepths = [];
    % Load heterogeneity config for Gauss functions that are called even if heterogeneity correction is turned off
    propHeterogeneity = matRad_HeterogeneityConfig();
end

% Check if correct base data is loaded for heterogeneity correction
if ~isempty(heteroCorrDepths) && ~isstruct(baseData.Z)
    matRad_cfg.dispWarning('calcParticleDoseBixel: heterogeneity correction enabled but no APM base data was loaded.')
end

% add potential offset
depths = baseData.depths + baseData.offset;

% convert from MeV cm^2/g per primary to Gy mm^2 per 1e6 primaries
conversionFactor = 1.6021766208e-02;


%% interpolate depth dose, sigmas and weights and calculate lateral sigmas
if isstruct(baseData.Z) && ~isfield(baseData,'sigma1')

    % interpolate depth dose and sigma
    X = matRad_interp1(depths,baseData.sigma,radDepths,'extrap');

    %compute lateral sigma
    sigmaSq = X.^2 + sigmaIni_sq;

elseif isstruct(baseData.Z) && isfield(baseData,'sigma1')

    % interpolate sigmas and weights
    X = matRad_interp1(depths,[baseData.sigma1 baseData.weight baseData.sigma2],radDepths);

    % compute lateral sigmas
    sigmaSq_Narr = X(:,1).^2 + sigmaIni_sq;
    sigmaSq_Bro  = X(:,3).^2 + sigmaIni_sq;

elseif ~isstruct(baseData.Z) && isfield(baseData,'sigma1')

    % interpolate depth dose, sigmas, and weights
    X = matRad_interp1(depths,[baseData.Z baseData.weight baseData.sigma1 baseData.sigma2],radDepths,'extrap');

    % set dose for query > tabulated depth dose values to zero
    X(radDepths > max(depths),1) = 0;

    % compute lateral sigmas
    sigmaSq_Narr = X(:,3).^2 + sigmaIni_sq;
    sigmaSq_Bro  = X(:,4).^2 + sigmaIni_sq;

else

    % interpolate depth dose and sigma
    X = matRad_interp1(depths,[baseData.Z baseData.sigma],radDepths);

    % set dose for query > tabulated depth dose values to zero
    X(radDepths > max(depths),1) = 0;

    %compute lateral sigma
    sigmaSq = X(:,2).^2 + sigmaIni_sq;
end

%  calculate lateral profiles
if isfield(baseData,'sigma1')

    L_Narr =  exp( -radialDist_sq ./ (2*sigmaSq_Narr))./(2*pi*sigmaSq_Narr); % Gauss
    L_Bro  =  exp( -radialDist_sq ./ (2*sigmaSq_Bro ))./(2*pi*sigmaSq_Bro );

    bixel.L = baseData.LatCutOff.CompFac * ((1-X(:,2)).*L_Narr + X(:,2).*L_Bro); % (1-w)*L_Narr + w*L_Bro

else

    bixel.L = baseData.LatCutOff.CompFac * exp( -radialDist_sq ./ (2*sigmaSq)) ./(2*pi*sigmaSq);

end

%% calculate sigma in range direction
if isstruct(baseData.Z)

    % calculate depthDoses with APM
    % no offset here...
    radDepths = radDepths - baseData.offset;

    % add sigma if heterogeneity correction wanted
    if ~isempty(heteroCorrDepths)
        switch propHeterogeneity.type
            case 'complete'
                [~,lungDepthAtBraggPeakIx] = min(abs(radialDist_sq+(radDepths-baseData.peakPos).^2));
                lungDepthAtBraggPeak = heteroCorrDepths(lungDepthAtBraggPeakIx);
                ellSq = ones(numel(radDepths),1)* (baseData.Z.width'.^2 + propHeterogeneity.getHeterogeneityCorrSigmaSq(lungDepthAtBraggPeak));

            case 'depthBased'
                % lungDepthAtGaussPeakIx = zeros(baseData.Z.mean)
                for i = 1:length(baseData.Z.mean)
                    [~,lungDepthAtGaussPeakIx(i)] = min(abs(radialDist_sq+(radDepths-baseData.Z.mean(i)).^2));
                end
                lungDepthAtGaussPeak = heteroCorrDepths(lungDepthAtGaussPeakIx);
                for i = 1:length(baseData.Z.mean)
                    ellSq(:,i) = ones(numel(radDepths),1)* (baseData.Z.width(i)'.^2 + propHeterogeneity.getHeterogeneityCorrSigmaSq(lungDepthAtGaussPeak(i)));
                end

            case 'voxelwise'
                ellSq = bsxfun(@plus, baseData.Z.width'.^2, propHeterogeneity.getHeterogeneityCorrSigmaSq(heteroCorrDepths, propHeterogeneity.modPower));
            %% Stolzenberg implementation:
            % local pmod for each bixel, mean of all pmod values for each
            % voxel in the bixel (see matRad_calcParticleDose.m)
            case 'local_pmod'
                if any(pmod~=0)
                    disp([heteroCorrDepths pmod])
                end
                ellSq = bsxfun(@plus, baseData.Z.width'.^2, propHeterogeneity.getHeterogeneityCorrSigmaSq(heteroCorrDepths, pmod));
            %%
            otherwise
                matRad_cfg.dispError('Error in heterogeneity correction')
        end
    else
        % use normal width to not include any depth modulation
        ellSq = ones(numel(radDepths),1)*baseData.Z.width'.^2;
    end

    % bixel.Z = (1./sqrt(2*pi*ellSq) .* exp(-bsxfun(@minus,baseData.Z.mean',radDepths).^2 ./ (2*ellSq)) )* baseData.Z.weight
    bixel.Z = propHeterogeneity.sumGauss(radDepths,baseData.Z.mean,ellSq',baseData.Z.weight);
else

    bixel.Z = X(:,1);

end

%% calculating the physical dose
bixel.physDose = conversionFactor * bixel.L .* bixel.Z;

%% manual modulation of alpha and sqrtBeta dose in case of RBE modeling
% Check if bioOptimization is needed (e.g. disabled for constRBE)
if propHeterogeneity.bioOpt
    % this only makes sense if APM data is used
    if isstruct(baseData.Z) && ~isempty(heteroCorrDepths)
        %% Stolzenberg implementation of modulation of alpha and sqrtBeta dose
        % This model is based on a linear Fit of initial beam energy (E), Depth
        % in Lung (WET) and modulation strength (Pmod)
        if propHeterogeneity.modulateRBE && strcmp(propHeterogeneity.type, 'local_pmod')
            
            % preallocate space for alpha beta
            tissueClasses = unique(vTissueIndex);
            if tissueClasses == 0
                matRad_cfg.dispError('Error in Heterogeneity correction: TissueIndex not assigned. ');
            end

            % Load fitparameters obtained from Monte Carlo simulations:
            load('Fit_alpha_beta_gaussian.mat');
            [A_alpha, mu_alpha, sigma_alpha] = propHeterogeneity.Biofitparam(fit.param_a,heteroCorrDepths./10,baseData.energy, pmod.*10^(-6));
            [A_beta, mu_beta, sigma_beta] = propHeterogeneity.Biofitparam(fit.param_b,heteroCorrDepths./10,baseData.energy, pmod.*10^(-6));
            
            bixel.Z_Aij = zeros(numel(radDepths),1);
            bixel.Z_Bij = zeros(numel(radDepths),1);

            % define gaussian:
            resolution   = min(diff(baseData.depths));
            bixel.heteroCorr.SigmaSq = propHeterogeneity.getHeterogeneityCorrSigmaSq(heteroCorrDepths,pmod);
            for i=1:length(radDepths)
                ix = vTissueIndex(i);  
                if bixel.heteroCorr.SigmaSq(i)~=0
                    %convolution needed: first define a finer dose grid to obtain the same resolution at each point 
                    bixel.heteroCorr.conv(i).fineGridradDepths = min(baseData.depths)-3*sqrt(bixel.heteroCorr.SigmaSq(i)):resolution:max(baseData.depths)+3*sqrt(bixel.heteroCorr.SigmaSq(i));
                    sigmaSq_eff = bsxfun(@plus, baseData.Z.width'.^2, bixel.heteroCorr.SigmaSq(i));
                    sigmaSq_eff = ones(numel(bixel.heteroCorr.conv(i).fineGridradDepths),1)*sigmaSq_eff;
                    bixel.heteroCorr.conv(i).Z_conv = propHeterogeneity.sumGauss(bixel.heteroCorr.conv(i).fineGridradDepths',baseData.Z.mean,sigmaSq_eff',baseData.Z.weight);
                    alpha_interp = matRad_interp1(baseData.depths, baseData.alpha, bixel.heteroCorr.conv(i).fineGridradDepths, 'extrap');
                    alphaDose = bixel.heteroCorr.conv(i).Z_conv.*alpha_interp;
                    beta_interp = matRad_interp1(baseData.depths, baseData.beta, bixel.heteroCorr.conv(i).fineGridradDepths, 'extrap');
                    sqrtbetaDose = bixel.heteroCorr.conv(i).Z_conv.*sqrt(beta_interp);
                    
                    % alpha convolution:
                    gaussNumOfPoints = round(6*sqrt(sigma_alpha(i))/resolution);
                    gaussWidth = linspace(-resolution*fix(gaussNumOfPoints/2), resolution*fix(gaussNumOfPoints/2), gaussNumOfPoints);
                    %save normalize Gaussian and grids
                    bixel.heteroCorr.conv(i).Gauss_alpha = propHeterogeneity.GaussBio(gaussWidth,A_alpha(i),mu_alpha(i),sigma_alpha(i));
                    % bixel.heteroCorr.conv(i).Gauss_alpha = bixel.heteroCorr.conv(i).Gauss_alpha/sum(bixel.heteroCorr.conv(i).Gauss_alpha);
                    
                    bixel.heteroCorr.conv(i).alphaDoseConv = conv(alphaDose(:,ix),bixel.heteroCorr.conv(i).Gauss_alpha,'same');
                    % Interpolation to radiological depth point:
                    bixel.Z_Aij(i) = conversionFactor * ...
                        matRad_interp1(bixel.heteroCorr.conv(i).fineGridradDepths,bixel.heteroCorr.conv(i).alphaDoseConv,radDepths(i),'extrap');
                    
                    % beta convolution:
                    gaussNumOfPoints = round(6*sqrt(sigma_beta(i))/resolution);
                    gaussWidth = linspace(-resolution*fix(gaussNumOfPoints/2), resolution*fix(gaussNumOfPoints/2), gaussNumOfPoints);
                    %save normalize Gaussian and grids
                    bixel.heteroCorr.conv(i).Gauss_beta = propHeterogeneity.GaussBio(gaussWidth,A_beta(i),mu_beta(i),sigma_beta(i));
                    % bixel.heteroCorr.conv(i).Gauss_beta = bixel.heteroCorr.conv(i).Gauss_beta/sum(bixel.heteroCorr.conv(i).Gauss_beta);
                    
                    bixel.heteroCorr.conv(i).betaDoseConv = conv(sqrtbetaDose(:,ix),bixel.heteroCorr.conv(i).Gauss_beta,'same');
                    % Interpolation to radiological depth point:
                    bixel.Z_Bij(i) = conversionFactor * ...
                        matRad_interp1(bixel.heteroCorr.conv(i).fineGridradDepths,bixel.heteroCorr.conv(i).betaDoseConv,radDepths(i),'extrap');
                    
                    
                else
                    % no convolution needed:
                    alphaDose = bixel.Z(i).*baseData.alpha(:,ix);
                    sqrtBetaDose = bixel.Z(i).*sqrt(baseData.beta(:,ix));

                    bixel.Z_Aij(i) = conversionFactor * ...
                        matRad_interp1(baseData.depths,alphaDose,radDepths(i),'extrap');
                    bixel.Z_Bij(i) = conversionFactor * ...
                        matRad_interp1(baseData.depths,sqrtBetaDose,radDepths(i),'extrap');
                end
            end
        else
            matRad_cfg.dispWarning('Modulation of RBE activated but no local modulation type is set.');
        end
        %% Homolka model:
        % Prepare heteroCorrSigma and Gaussian convolution
        if propHeterogeneity.modulateLET || propHeterogeneity.modulateBioDose

            [~,lungDepthAtBraggPeakIx] = min(abs(radialDist_sq+(radDepths-baseData.peakPos).^2));
            lungDepthAtBraggPeak = heteroCorrDepths(lungDepthAtBraggPeakIx);
            bixel.heteroCorr.SigmaSq = propHeterogeneity.getHeterogeneityCorrSigmaSq(lungDepthAtBraggPeak);

            % define Gaussian around zero
            resolution   = min(diff(baseData.depths));
            gaussNumOfPoints = round(6*sqrt(bixel.heteroCorr.SigmaSq)/resolution);
            gaussWidth = linspace(-resolution*fix(gaussNumOfPoints/2), resolution*fix(gaussNumOfPoints/2), gaussNumOfPoints);
            %save normalize Gaussian and grids
            bixel.heteroCorr.Gauss = propHeterogeneity.Gauss(gaussWidth,0,bixel.heteroCorr.SigmaSq);
            bixel.heteroCorr.Gauss = bixel.heteroCorr.Gauss/sum(bixel.heteroCorr.Gauss);
            bixel.heteroCorr.fineGrid = min(baseData.depths)-3*sqrt(bixel.heteroCorr.SigmaSq):resolution:max(baseData.depths)+3*sqrt(bixel.heteroCorr.SigmaSq);
            bixel.heteroCorr.coarseGrid = baseData.depths;

        end

        % LET convolution
        if propHeterogeneity.modulateLET && isfield(baseData,'LET')

            % LET extrapolation to finer grid
            LET = matRad_interp1(bixel.heteroCorr.coarseGrid,baseData.LET,bixel.heteroCorr.fineGrid,'extrap');
            LET(isnan(LET)) = 0;

            if bixel.heteroCorr.SigmaSq ~= 0
                % Convolution
                bixel.LET = conv(LET,bixel.heteroCorr.Gauss,'same');

                % get values for individual radDepths
                bixel.LET = matRad_interp1(bixel.heteroCorr.fineGrid,bixel.LET,radDepths,'extrap');
            else
                bixel.LET = matRad_interp1(bixel.heteroCorr.fineGrid,LET,radDepths,'extrap');
            end
            %      figure, plot(x,LET), hold on, plot(x,bixel.LET)
        end

        if propHeterogeneity.modulateBioDose
            % preallocate space for alpha beta
            tissueClasses = unique(vTissueIndex);
            if tissueClasses == 0
                matRad_cfg.dispError('Error in Heterogeneity correction: TissueIndex not assigned. ');
            end

            bixel.Z_Aij = zeros(numel(radDepths),1);
            bixel.Z_Bij = zeros(numel(radDepths),1);

            if isfield(baseData,'alphaDose') && isstruct(baseData.alphaDose)
                for i = 1:numel(tissueClasses)
                    ix = vTissueIndex == tissueClasses(i);
                    bixel.Z_Aij(ix)  = conversionFactor * ...
                        propHeterogeneity.sumGauss(radDepths(ix),baseData.alphaDose(tissueClasses(i)).mean,...
                        (baseData.alphaDose(tissueClasses(i)).width).^2 + bixel.heteroCorr.SigmaSq, ...
                        baseData.alphaDose(tissueClasses(i)).weight);

                    bixel.Z_Bij(ix)  = conversionFactor * ...
                        propHeterogeneity.sumGauss(radDepths(ix),baseData.SqrtBetaDose(tissueClasses(i)).mean,...
                        (baseData.SqrtBetaDose(tissueClasses(i)).width).^2 + bixel.heteroCorr.SigmaSq, ...
                        baseData.SqrtBetaDose(tissueClasses(i)).weight)';
                end

            elseif isfield(baseData,'alpha')
                alphaDose = baseData.Z.profileORG .* baseData.alpha;
                sqrtBetaDose = baseData.Z.profileORG .* sqrt(baseData.beta);

                % alpha beta convolution
                alphaDose = matRad_interp1(bixel.heteroCorr.coarseGrid,alphaDose,bixel.heteroCorr.fineGrid,'extrap');
                alphaDose(isnan(alphaDose)) = 0;
                sqrtBetaDose = matRad_interp1(bixel.heteroCorr.coarseGrid,sqrtBetaDose,bixel.heteroCorr.fineGrid,'extrap');
                sqrtBetaDose(isnan(sqrtBetaDose)) = 0;

                if bixel.heteroCorr.SigmaSq ~= 0
                    alphaDoseConv = zeros(size(alphaDose));
                    sqrtBetaDoseConv = zeros(size(alphaDose));
                    for i = 1:size(alphaDose,2)
                        alphaDoseConv(:,i) = conv(alphaDose(:,i),bixel.heteroCorr.Gauss,'same');
                        sqrtBetaDoseConv(:,i) = conv(sqrtBetaDose(:,i),bixel.heteroCorr.Gauss,'same');
                    end
                else
                    alphaDoseConv = alphaDose;
                    sqrtBetaDoseConv = sqrtBetaDose;
                end

                for i = 1:numel(tissueClasses)
                    ix = vTissueIndex == tissueClasses(i);
                    bixel.Z_Aij(ix) = conversionFactor * ...
                        matRad_interp1(bixel.heteroCorr.fineGrid,alphaDoseConv(:,tissueClasses(i)),radDepths(ix),'extrap');
                    bixel.Z_Bij(ix) = conversionFactor * ...
                        matRad_interp1(bixel.heteroCorr.fineGrid,sqrtBetaDoseConv(:,tissueClasses(i)),radDepths(ix),'extrap');
                end

            else
                matRad_cfg.dispError('Error in RBE bixel calculation.');
            end

            if any(isnan(bixel.Z_Aij(:))) || all(bixel.Z_Aij(:) == 0)
                matRad_cfg.dispError('unexpected NaN or all zero detected in heterogeneity convolution.')
            end
        end

    end
end

% check if we have valid dose values
if any(isnan(bixel.physDose)) || any(bixel.physDose<0)
    error('Error in particle dose calculation.');
end
