
matRad_cfg =  MatRad_Config.instance();

% default: dose influence matrix computation
if ~exist('calcDoseDirect','var')
    calcDoseDirect = false;
end

% to guarantee downwards compatibility with data that does not have
% ct.x/y/z
if ~any(isfield(ct,{'x','y','z'}))
    ct.x = ct.resolution.x*[0:ct.cubeDim(2)-1]-ct.resolution.x/2;
    ct.y = ct.resolution.y*[0:ct.cubeDim(1)-1]-ct.resolution.y/2;
    ct.z = ct.resolution.z*[0:ct.cubeDim(3)-1]-ct.resolution.z/2;
end

% set grids
if ~isfield(pln,'propDoseCalc') || ...
   ~isfield(pln.propDoseCalc,'doseGrid') || ...
   ~isfield(pln.propDoseCalc.doseGrid,'resolution')
    % default values
    dij.doseGrid.resolution = matRad_cfg.propDoseCalc.defaultResolution;
else
    % take values from pln strcut
    dij.doseGrid.resolution.x = pln.propDoseCalc.doseGrid.resolution.x;
    dij.doseGrid.resolution.y = pln.propDoseCalc.doseGrid.resolution.y;
    dij.doseGrid.resolution.z = pln.propDoseCalc.doseGrid.resolution.z;
end

dij.doseGrid.x = ct.x(1):dij.doseGrid.resolution.x:ct.x(end);
dij.doseGrid.y = ct.y(1):dij.doseGrid.resolution.y:ct.y(end);
dij.doseGrid.z = ct.z(1):dij.doseGrid.resolution.z:ct.z(end);

dij.doseGrid.dimensions  = [numel(dij.doseGrid.y) numel(dij.doseGrid.x) numel(dij.doseGrid.z)];
dij.doseGrid.numOfVoxels = prod(dij.doseGrid.dimensions);

dij.ctGrid.resolution.x = ct.resolution.x;
dij.ctGrid.resolution.y = ct.resolution.y;
dij.ctGrid.resolution.z = ct.resolution.z;

dij.ctGrid.x = ct.x;
dij.ctGrid.y = ct.y;
dij.ctGrid.z = ct.z;

dij.ctGrid.dimensions  = [numel(dij.ctGrid.y) numel(dij.ctGrid.x) numel(dij.ctGrid.z)];
dij.ctGrid.numOfVoxels = prod(dij.ctGrid.dimensions);

% adjust isocenter internally for different dose grid
offset = [dij.doseGrid.resolution.x - dij.ctGrid.resolution.x ...
          dij.doseGrid.resolution.y - dij.ctGrid.resolution.y ...
          dij.doseGrid.resolution.z - dij.ctGrid.resolution.z];
    
for i = 1:numel(stf)
    stf(i).isoCenter = stf(i).isoCenter + offset;
end

%set up HU to rED or rSP conversion
if ~isfield(pln,'propDoseCalc') || ~isfield(pln.propDoseCalc,'useGivenEqDensityCube')
    disableHUconversion = matRad_cfg.propDoseCalc.defaultUseGivenEqDensityCube;
else
    disableHUconversion = pln.propDoseCalc.useGivenEqDensityCube;
end

%If we want to omit HU conversion check if we have a ct.cube ready
if disableHUconversion && ~isfield(ct,'cube')
    matRad_cfg.dispWarning('HU Conversion requested to be omitted but no ct.cube exists! Will override and do the conversion anyway!');
    disableHUconversion = false;
end
    
% calculate rED or rSP from HU
if disableHUconversion
    matRad_cfg.dispInfo('Omitting HU to rED/rSP conversion and using existing ct.cube!\n');
else
    ct = matRad_calcWaterEqD(ct, pln);
end

% meta information for dij
dij.numOfBeams         = pln.propStf.numOfBeams;
dij.numOfScenarios     = pln.multScen.numOfCtScen;
dij.numOfRaysPerBeam   = [stf(:).numOfRays];
dij.totalNumOfBixels   = sum([stf(:).totalNumOfBixels]);
dij.totalNumOfRays     = sum(dij.numOfRaysPerBeam);

% check if full dose influence data is required
if calcDoseDirect 
    numOfColumnsDij      = length(stf);
    numOfBixelsContainer = 1;
else
    numOfColumnsDij      = dij.totalNumOfBixels;
    numOfBixelsContainer = ceil(dij.totalNumOfBixels/10);
end

% set up arrays for book keeping
dij.bixelNum = NaN*ones(numOfColumnsDij,1);
dij.rayNum   = NaN*ones(numOfColumnsDij,1);
dij.beamNum  = NaN*ones(numOfColumnsDij,1);


% Allocate space for dij.physicalDose sparse matrix
for ctScen = 1:pln.multScen.numOfCtScen
   for shiftScen = 1:pln.multScen.totNumShiftScen
      for rangeShiftScen = 1:pln.multScen.totNumRangeScen   
         if pln.multScen.scenMask(ctScen,shiftScen,rangeShiftScen)
            dij.physicalDose{ctScen,shiftScen,rangeShiftScen} = spalloc(dij.doseGrid.numOfVoxels,numOfColumnsDij,1);
         end
      end
   end
end

% Allocate memory for dose_temp cell array
doseTmpContainer = cell(numOfBixelsContainer,pln.multScen.numOfCtScen,pln.multScen.totNumShiftScen,pln.multScen.totNumRangeScen);

% take only voxels inside patient or as specified in
% pln.propDoseCalc.voxelSubIx
if ~isfield(pln,'propDoseCalc') || ~isfield(pln.propDoseCalc,'voxelSubIx')
    subIx = matRad_cfg.propDoseCalc.defaultVoxelSubIx;
else
    subIx = pln.propDoseCalc.voxelSubIx;
end

if isempty(subIx)
    VctGrid = [cst{:,4}];
    VctGrid = unique(vertcat(VctGrid{:}));
else
    VctGrid = subIx;
end

% ignore densities outside of contours
if ~isfield(pln,'propDoseCalc') || ~isfield(pln.propDoseCalc,'ignoreOutsideDensities')
    ignoreOutsideDensities = matRad_cfg.propDoseCalc.defaultIgnoreOutsideDensities;
else
    ignoreOutsideDensities = pln.propDoseCalc.ignoreOutsideDensities;
end

if ignoreOutsideDensities
    eraseCtDensMask = ones(prod(ct.cubeDim),1);
    eraseCtDensMask(VctGrid) = 0;
    for i = 1:ct.numOfCtScen
        ct.cube{i}(eraseCtDensMask == 1) = 0;
    end
end

% ser overlap prioriites
cst = matRad_setOverlapPriorities(cst);

% resizing cst to dose cube resolution
cst = matRad_resizeCstToGrid(cst,dij.ctGrid.x,dij.ctGrid.y,dij.ctGrid.z,...
   dij.doseGrid.x,dij.doseGrid.y,dij.doseGrid.z);

% Convert CT subscripts to linear indices.
[yCoordsV_vox, xCoordsV_vox, zCoordsV_vox] = ind2sub(ct.cubeDim,VctGrid);

% receive linear indices and grid locations from the dose grid
tmpCube    = zeros(ct.cubeDim);
tmpCube(VctGrid) = 1;
% interpolate cube
VdoseGrid = find(matRad_interp3(dij.ctGrid.x,  dij.ctGrid.y,   dij.ctGrid.z,tmpCube, ...
                                dij.doseGrid.x,dij.doseGrid.y',dij.doseGrid.z,'nearest'));

% Convert CT subscripts to coarse linear indices.
[yCoordsV_voxDoseGrid, xCoordsV_voxDoseGrid, zCoordsV_voxDoseGrid] = ind2sub(dij.doseGrid.dimensions,VdoseGrid);

% load base data
fileName = [pln.radiationMode '_' pln.machine];
try
   load([fileparts(mfilename('fullpath')) filesep 'basedata' filesep fileName '.mat']);
catch
   matRad_cfg.dispError('Could not find the following machine file: %s\n',fileName); 
end

% compute SSDs -> Removed for now because it is scenario-dependent
% stf = matRad_computeSSD(stf,ct);

%% lungmodulation implementation
%  create a variable analogue to "ct" for "lung" to hand to raytracer
if isfield(pln, 'lungModulation') && pln.lungModulation.CalcHetero
    % modulation cube with properties of ct      
    modulation = ct;
    modulation.numOfCtScen = 2;
    % cube{1} = logical Mask of all Voxels inside Lung / or rather
    % modulating voxels (see HU trehshold later)
    % cube{2} = all voxels inside lung with its final SPR assigned all
    % other voxels = 0
    modulation.cube(2) = {zeros(ct.cubeDim)};
    % cube{1} = logical mask with 1 inside lung
    modulation.cube(1) = {zeros(ct.cubeDim)};
    modulation.cubeHU = {ones(ct.cubeDim).*-1024};
    % in the future: set additional filter for modulating material atm its just "lung"
    for innerloop = 1:size(cst,1)
        %             if strcmpi(cst{innerloop,2}, 'lung')
        if contains(cst{innerloop,2}, 'Lung') || contains(cst{innerloop,2}, 'lung')

            modulation.cube{2}(cst{innerloop,4}{1}) = ct.cube{1}(cst{innerloop,4}{1});
            modulation.cubeHU{1}(cst{innerloop,4}{1}) = ct.cubeHU{1}(cst{innerloop,4}{1});
            %                 % take only voxels inside lung
            %                 VlungGrid = [cst{innerloop,4}];
            %                 VlungGrid = unique(vertcat(VlungGrid{:}));
        end
%         assignin('base', 'modulation', modulation)
    end
    % apply HU threshold 
    HU_schwelle = [-900 -100];
    modulation.cube{2}(modulation.cubeHU{1}<HU_schwelle(1) | modulation.cubeHU{1}>HU_schwelle(2)) = 0;
    modulation.cube{1}(modulation.cube{2}>0) = 1;
    assignin('base', 'modulation', modulation)


%     filterkernel_std = ones(3,3,3);
%     filterkernel_mean = ones(3,3,3)./(3^3);
%     standardDevImage = stdfilt(modulation.cubeHU{1},filterkernel_std);
%     meanImage = convn(modulation.cubeHU{1}, filterkernel_mean, 'same');
    
%     % EQ7 Phys. Med. Biol. 66 (2021) 185002
%     Pfit = (standardDevImage.^2 ./ (-1000 .* meanImage - meanImage.^2) );
%     modulation.cube(1) = {nthroot( Pfit,3)};
%     modulation.cube{1}(modulation.cube{1}<0) = 0;
%     modulation.cube{1}(modulation.cube{1}>1) = 1;
% %     temporary mask with 1 inside lung
%     modulation.cube{1}(modulation.cube{2}>0) = 1;
% %  additional filter for HU range
%     modulation.cube{1}(modulation.cubeHU{1}<-900 & modulation.cubeHU{1}>-100) = 0)

    figure,imshow(modulation.cube{2}(:,:,round(modulation.cubeDim(3)/2)))
end