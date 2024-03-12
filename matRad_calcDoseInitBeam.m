matRad_cfg.dispInfo('Beam %d of %d:\n',i,dij.numOfBeams);

% remember beam and bixel number
if calcDoseDirect
    dij.beamNum(i)    = i;
    dij.rayNum(i)     = i;
    dij.bixelNum(i)   = i;
end

bixelsPerBeam = 0;

% convert voxel indices to real coordinates using iso center of beam i
xCoordsV       = xCoordsV_vox(:)*ct.resolution.x-stf(i).isoCenter(1);
yCoordsV       = yCoordsV_vox(:)*ct.resolution.y-stf(i).isoCenter(2);
zCoordsV       = zCoordsV_vox(:)*ct.resolution.z-stf(i).isoCenter(3);
coordsV        = [xCoordsV yCoordsV zCoordsV];

xCoordsVdoseGrid = xCoordsV_voxDoseGrid(:)*dij.doseGrid.resolution.x-stf(i).isoCenter(1);
yCoordsVdoseGrid = yCoordsV_voxDoseGrid(:)*dij.doseGrid.resolution.y-stf(i).isoCenter(2);
zCoordsVdoseGrid = zCoordsV_voxDoseGrid(:)*dij.doseGrid.resolution.z-stf(i).isoCenter(3);
coordsVdoseGrid  = [xCoordsVdoseGrid yCoordsVdoseGrid zCoordsVdoseGrid];

% Get Rotation Matrix
% Do not transpose matrix since we usage of row vectors &
% transformation of the coordinate system need double transpose

rotMat_system_T = matRad_getRotationMatrix(stf(i).gantryAngle,stf(i).couchAngle);

% Rotate coordinates (1st couch around Y axis, 2nd gantry movement)
rot_coordsV         = coordsV*rotMat_system_T;
rot_coordsVdoseGrid = coordsVdoseGrid*rotMat_system_T;

rot_coordsV(:,1) = rot_coordsV(:,1)-stf(i).sourcePoint_bev(1);
rot_coordsV(:,2) = rot_coordsV(:,2)-stf(i).sourcePoint_bev(2);
rot_coordsV(:,3) = rot_coordsV(:,3)-stf(i).sourcePoint_bev(3);

rot_coordsVdoseGrid(:,1) = rot_coordsVdoseGrid(:,1)-stf(i).sourcePoint_bev(1);
rot_coordsVdoseGrid(:,2) = rot_coordsVdoseGrid(:,2)-stf(i).sourcePoint_bev(2);
rot_coordsVdoseGrid(:,3) = rot_coordsVdoseGrid(:,3)-stf(i).sourcePoint_bev(3);

% calculate geometric distances
geoDistVdoseGrid{1}= sqrt(sum(rot_coordsVdoseGrid.^2,2));

% Calculate radiological depth cube
matRad_cfg.dispInfo('matRad: calculate radiological depth cube... ');
if strcmp(pln.propDoseCalc.fineSampling.calcMode, 'fineSampling')
    [radDepthVctGrid, radDepthsMat] = matRad_rayTracing(stf(i),ct,VctGrid,rot_coordsV,pln.propDoseCalc.effectiveLateralCutOff);
else
    radDepthVctGrid = matRad_rayTracing(stf(i),ct,VctGrid,rot_coordsV,pln.propDoseCalc.effectiveLateralCutOff);
end
matRad_cfg.dispInfo('done.\n');

% interpolate radiological depth cube to dose grid resolution
radDepthVdoseGrid = matRad_interpRadDepth(ct,VctGrid,VdoseGrid,dij.doseGrid.x,dij.doseGrid.y,dij.doseGrid.z,radDepthVctGrid);
 
if exist('radDepthsMat', 'var')
    % interpolate radiological depth cube used for fine sampling to dose grid resolution    
    radDepthsMat = cellfun(@(radDepthCube) matRad_interp3(dij.ctGrid.x,  dij.ctGrid.y,   dij.ctGrid.z,radDepthCube,dij.doseGrid.x,dij.doseGrid.y',dij.doseGrid.z,'nearest'),radDepthsMat,'UniformOutput',false);    
end

% limit rotated coordinates to positions where ray tracing is availabe
rot_coordsVdoseGrid = rot_coordsVdoseGrid(~isnan(radDepthVdoseGrid{1}),:);

%% lungmodulation implementation
% calculate raytracing for all voxels inside lung
if isfield(pln, 'lungModulation') && pln.lungModulation.CalcHetero
    matRad_cfg.dispInfo('matRad: calculate depth cube for lungMod...');
    modulationDepthVctGrid = matRad_rayTracing(stf(i),modulation,VctGrid,rot_coordsV,pln.propDoseCalc.effectiveLateralCutOff);
%     lungDepthVctGrid = matRad_rayTracing(stf(i),modulation,VlungGrid,rot_coordsV,effectiveLateralCutoff);
    matRad_cfg.dispInfo('done.\n');
    % interpolate modulation depth cube to dose grid resolution
    
    modulationDepthVdoseGrid = matRad_interpRadDepth...
    (modulation,VctGrid,VdoseGrid,dij.doseGrid.x,dij.doseGrid.y,dij.doseGrid.z,modulationDepthVctGrid);
%         lungDepthVdoseGrid(currentScen) = matRad_interpRadDepth...
%             (modulation,currentScen,VctGrid,VdoseGrid,dij.doseGrid.x,dij.doseGrid.y,dij.doseGrid.z,modulationDepthVctGrid);
     
    if ~isfield(pln.lungModulation, 'ModPower')
        pln.lungModulation.ModPower = 250;
        matRad_cfg.dispInfo('matRad: no modulation power specified using default value of Pmod = 250\mum ...');
    end
%   prepare basedata for the convolution => 
%   precomputation of prolonged basedata needed for
%   convolution in calcParticleDoseBixel
%   currently all energies are prolonged => consider prolonging only used
%   energies
% min_energyIx = min(round2(stf(i).ray.energy,4)) == round2([machine.data.energy],4);
% min_energyIx = find(min(round2(stf.ray.energy,4)) == round2([machine.data.energy],4));
min_energyIx = 1;
% max_energyIx = max(round2(stf(i).ray(j).energy,4)) == round2([machine.data.energy],4);
% max_energyIx  = find(max(round2(stf.ray.energy,4)) == round2([machine.data.energy],4));
max_energyIx = size(machine.data,2);

for energiestep = min_energyIx:max_energyIx
    [machine.data(energiestep).Z_adapted, machine.data(energiestep).depths_adapted] = matRad_extendBaseData(...
        machine.data(energiestep));
end

end
 