function vol = ldrb_main_adapted(config)
% Laplace–Dirichlet Rule-Based (LDRB) algorithm for ventricular fibers
% from Bayer 2012: https://doi.org/10.1007/s10439-012-0593-5
% (see algorithm "DefineFibers" in the supplement),
%
% The original algorithm was adapted to yield a continuous and
% approximately linear rotation of fibers across the free walls.

% This function can be applied to BI- and ONE-ventricular geometries.

if ~isfield(config,'sourceDir')
    error('Mandatory parameter ''sourceDir'' missing in config.');
end
if ((isfield(config,'exportFinalResult') && config.exportFinalResult) || ...
    (isfield(config,'exportIntermediateResults') && config.exportIntermediateResults)) && ...
   ~isfield(config,'targetPrefix')
    error('Mandatory parameter ''targetPrefix'' missing in config.');
end

% default parameters
cfg.sourceDir = '';
cfg.targetPrefix = '';
cfg.onlyOneVentricle = false;
cfg.volName = 'heart';
cfg.surNames = {'epi','lv','rv','base','apex'};
cfg.alphaSeptLeft  = 60;
cfg.alphaSeptRight = 60;
cfg.alphaWallEndo  = 60;
cfg.alphaWallEpi   = -60;
cfg.betaSeptLeft   = -45;
cfg.betaSeptRight  = -45;
cfg.betaWallEndo   = -45;
cfg.betaWallEpi    = 45;
cfg.exportIntermediateResults = false;
cfg.exportFinalResult = true;
cfg.exportFiber = true;
cfg.exportSheet = true;
cfg.exportSheetnormal = true;
cfg.exportAngles = true;
cfg.outputAngleUnit = 'rad'; % 'rad', 'deg' or 'ibt'
cfg.exportDebugAngle = false;
cfg.tol = 1e-12;
cfg.maxit = 1000;

% overwrite default parameters, if provided by the user
fns = fieldnames(config);
for i = 1:numel(fns)
    fn = fns{i};
    if ~isfield(cfg,fn)
        error('Unknown parameter ''%s'' found in config.', fn);
    end
    cfg.(fn) = config.(fn);
end

%% Convert angles to radian

alphaSeptLeft  = deg2rad(cfg.alphaSeptLeft);
alphaSeptRight = deg2rad(cfg.alphaSeptRight);
alphaWallEndo  = deg2rad(cfg.alphaWallEndo);
alphaWallEpi   = deg2rad(cfg.alphaWallEpi);
betaSeptLeft   = deg2rad(cfg.betaSeptLeft);
betaSeptRight  = deg2rad(cfg.betaSeptRight);
betaWallEndo   = deg2rad(cfg.betaWallEndo);
betaWallEpi    = deg2rad(cfg.betaWallEpi);

%%
if cfg.onlyOneVentricle
    fprintf('\n==== Generating fibers for one-ventricular geometry ====\n\n');
else
    fprintf('\n==== Generating fibers for bi-ventricular geometry ====\n\n');
end

%% Check if target file can be written

if cfg.exportFinalResult || cfg.exportIntermediateResults
    targetFile = sprintf('%s.vtu', cfg.targetPrefix);
    if exist(targetFile, 'file')
        fid = fopen(targetFile, 'r+');
    else
        fid = fopen(targetFile, 'w');
    end
	if fid == -1
		error('Could not write to ''%s.vtu''.', cfg.targetPrefix);
	end
	fclose(fid);
end

%% Set up boundary conditions for Laplace equation from surfaces

fprintf('Setting up boundary conditions...           '); tic;

vol = vtkRead(sprintf('%s/%s.vtu', cfg.sourceDir, cfg.volName));
res = vol;
numCells = size(vol.cells,1);

ids = cell(5,1);
if numel(cfg.surNames) == 1
    sur = vtkRead(sprintf('%s/%s.vtp', cfg.sourceDir, cfg.surNames{1}));
    surToVol = vtkMapPointIds(vol, sur);
    ids{1} = unique(surToVol(sur.pointData.class==2 | sur.pointData.class==5)); % epi
    ids{2} = unique(surToVol(sur.pointData.class==3)); % lv
    ids{3} = unique(surToVol(sur.pointData.class==4)); % rv
    ids{4} = unique(surToVol(sur.pointData.class==1)); % base
    ids{5} = unique(surToVol(sur.pointData.class==5)); % apex
elseif numel(cfg.surNames) >= 4
    for i = 1:numel(cfg.surNames)
        list = dir(sprintf('%s/%s.*', cfg.sourceDir, cfg.surNames{i}));
        sur = [];
        for j = 1:numel(list)
            name = list(j).name;
            if endsWith(name, {'.stl','.ply','.obj','.vtp'})
                sur = vtkRead(sprintf('%s/%s', cfg.sourceDir, name));
                break;
            end
        end
        if isempty(sur)
            error('Could not find a surface file with the name ''%s'' in the folder ''%s''.', cfg.surNames{i}, cfg.sourceDir);
        end
        ids{i} = unique(vtkMapPointIds(vol, sur));
    end
else
    error('cfg.surNames must contain either 1 name (vtp with point data ''class'') or >=4 names (individual surface files for each region).')
end

if cfg.onlyOneVentricle
    % make sure pairs of ids are unique
    ids{1} = setdiff(ids{1}, ids{2});
    ids{3} = setdiff(ids{3}, ids{4});
    % lv boundary conditions
    ids_lv = [ids{1}; ids{2}];
    val_lv = [zeros(size(ids{1})); ones(size(ids{2}))];
    % apico-basal boundary conditions
    ids_ab = [ids{3}; ids{4}];
    val_ab = [ones(size(ids{3})); zeros(size(ids{4}))];
else
    % make sure pairs of ids are unique
    ids{1} = setdiff(ids{1}, [ids{2}; ids{3}]);
    ids{2} = setdiff(ids{2}, ids{3});
    ids{4} = setdiff(ids{4}, ids{5});
    % lv boundary conditions
    ids_lv = [ids{1}; ids{2}; ids{3}];
    val_lv = [zeros(size(ids{1})); ones(size(ids{2})); zeros(size(ids{3}))];
    % rv boundary conditions
    ids_rv = [ids{1}; ids{2}; ids{3}];
    val_rv = [zeros(size(ids{1})); zeros(size(ids{2})); ones(size(ids{3}))];
    % apico-basal boundary conditions
    ids_ab = [ids{4}; ids{5}];
    val_ab = [ones(size(ids{4})); zeros(size(ids{5}))];
end

fprintf('%.1f seconds\n', toc);

%% Computing differential operators

fprintf('Computing differential operators...         '); tic;

P = double(vol.points); % points
C = double(vol.cells);  % cells
L = cotmatrix(P,C);     % Laplacian
G = grad(P,C);          % gradient

fprintf('%.1f seconds\n', toc);

%% Compute Laplace solutions

fprintf('Computing Laplace solutions...              '); tic

lapLv = solveLaplace(L, ids_lv, val_lv, cfg.tol, cfg.maxit);
lapAb = solveLaplace(L, ids_ab, val_ab, cfg.tol, cfg.maxit);

res.pointData.lapLv = lapLv;
res.pointData.lapAb = lapAb;

if ~cfg.onlyOneVentricle
   lapRv = solveLaplace(L, ids_rv, val_rv, cfg.tol, cfg.maxit);
   lapEpi = 1 - (lapLv + lapRv);
   
   res.pointData.lapRv = lapRv;
   res.pointData.lapEpi = lapEpi;
end

% vol.pointData.laplace_ab = single(lapAb);
% vol.pointData.laplace_epi = single(lapEpi);
% vol.pointData.laplace_lv = single(lapLv);
% vol.pointData.laplace_rv = single(lapRv);

fprintf('%.1f seconds\n', toc);

%% Compute gradients of Laplace solutions

fprintf('Computing gradients of Laplace solutions... '); tic;

if ~cfg.onlyOneVentricle
    gradEpi = reshape(G*lapEpi, size(C,1), 3);
    gradRv  = reshape(G*lapRv, size(C,1), 3); 
    gradEpi = normalizeRows(gradEpi, cfg.tol);
    gradRv  = normalizeRows(gradRv, cfg.tol);
    
    res.cellData.gradEpi = gradEpi;
    res.cellData.gradRv = gradRv;
else
    gradEpi = NaN(numCells,1);
    gradRv = NaN(numCells,1);
end

gradLv = reshape(G*lapLv, size(C,1), 3);
gradAb = reshape(G*lapAb, size(C,1), 3);
gradLv = normalizeRows(gradLv, cfg.tol);
gradAb = normalizeRows(gradAb, cfg.tol);

res.cellData.gradLv = gradLv;
res.cellData.gradAb = gradAb;

fprintf('%.1f seconds\n', toc);

%% Convert Laplace solutions from point to cell data

if ~cfg.onlyOneVentricle
    lapEpi = mean(lapEpi(C),2);
    lapRv  = mean(lapRv(C),2);
    
    res.cellData.lapEpi = lapEpi;
    res.cellData.lapRv = lapRv;
end

lapLv = mean(lapLv(C),2);
lapAb = mean(lapAb(C),2);

res.cellData.lapLv = lapLv;
res.cellData.lapAb = lapAb;

%% Export intermediate results

if cfg.exportIntermediateResults
    fprintf('Exporting intermediate results...           '); tic;

    vtkWrite(res, sprintf('%s_intermediateResults.vtu', cfg.targetPrefix));

    fprintf('%.1f seconds\n', toc);
end

%% Compute local angles

if cfg.onlyOneVentricle
    tLeftRight = NaN(numCells,1);
    alphaSept = NaN(numCells,1);
    betaSept = NaN(numCells,1);
    
    tEndoEpi = lapLv;
    alphaWall = (1-tEndoEpi) * alphaWallEndo + tEndoEpi * alphaWallEpi;
    betaWall  = (1-tEndoEpi) * betaWallEndo  + tEndoEpi * betaWallEpi;
else
    tLeftRight = lapRv./(lapLv+lapRv);
    tLeftRight(lapLv+lapRv < cfg.tol) = 0.5;
    alphaSept = (1-tLeftRight) * alphaSeptLeft - tLeftRight * alphaSeptRight;
    betaSept  = (1-tLeftRight) * betaSeptLeft  - tLeftRight * betaSeptRight;
    
    tEndoEpi = sqrt(max(lapEpi,0)); % adapted compared to the original algorithm
    alphaWall = (1-tEndoEpi) * alphaWallEndo + tEndoEpi * alphaWallEpi;
    betaWall  = (1-tEndoEpi) * betaWallEndo  + tEndoEpi * betaWallEpi;
end

%% Compute fibers for each cell

fprintf('Computing fibers for each cell...                       '); tic;

fiber = NaN(numCells,3,'single');
sheet = fiber;
sheetnormal = fiber;
theta = NaN(numCells,1,'single');
phi   = theta;

dq = parallel.pool.DataQueue;
progressStep = floor(numCells/100);
progressCount = 0;
afterEach(dq, @printProgress);
function printProgress(~)
    progressCount = progressCount+1;
    fprintf('%s%03i percent ', repmat(char(8),12,1), 10*round(progressCount/10));
end

parfor i = 1:numCells
    
    if cfg.onlyOneVentricle
        R     = ldrb_orient(ldrb_axis(gradAb(i,:)', gradLv(i,:)', cfg.tol), alphaWall(i), betaWall(i));
    else
        Qlv   = ldrb_orient(ldrb_axis(gradAb(i,:)', -gradLv(i,:)', cfg.tol), alphaSept(i), betaSept(i));
        Qrv   = ldrb_orient(ldrb_axis(gradAb(i,:)',  gradRv(i,:)', cfg.tol), alphaSept(i), betaSept(i));
        Qendo = ldrb_bislerp_adapted(Qlv, Qrv, tLeftRight(i)); % adapted compared to the original algorithm
        
        % adapted compared to the original algorithm
        if tLeftRight(i) > 0.5
            Qendo(:,1) = -Qendo(:,1);
            Qendo(:,3) = -Qendo(:,3);
        end
        
        Qepi  = ldrb_orient(ldrb_axis(gradAb(i,:)', gradEpi(i,:)', cfg.tol), alphaWall(i), betaWall(i));
        R     = ldrb_bislerp_adapted(Qendo, Qepi, tEndoEpi(i)); % adapted compared to the original algorithm
    end
    
    f  = R(:,1)';
    sn = R(:,2)';
    s  = R(:,3)';
    fiber(i,:) = f;
    sheet(i,:) = s;
    sheetnormal(i,:) = sn;
    
    % compute fiber angles
    t = acos(f(3));
    p = atan2(f(2), f(1));
    if p < 0
        t = acos(-f(3));
        p = atan2(-f(2), -f(1));
    end
    switch cfg.outputAngleUnit
        case 'deg' 
            t = t*180/pi;
            p = p*180/pi;
        case 'ibt'
            t = round(t*254/pi);
            p = round(p*254/pi);
    end
    theta(i) = t;
    phi(i) = p;
    
    if ~mod(i,progressStep)
        send(dq, 0);
    end

end

if strcmp(cfg.outputAngleUnit, 'ibt')
    theta = uint8(theta);
    phi = uint8(phi);
end

res.cellData.Fiber = fiber;
res.cellData.Sheet = sheet;
res.cellData.Sheetnormal = sheetnormal;
res.cellData.Theta = theta;
res.cellData.Phi = phi;

fprintf('%s%.1f seconds\n', repmat(char(8),12,1), toc);

%% Exporting result

if cfg.exportFiber
    vol.cellData.Fiber = fiber;
end
if cfg.exportSheet
    vol.cellData.Sheet = sheet;
end
if cfg.exportSheetnormal
    vol.cellData.Sheetnormal = sheetnormal;
end
if cfg.exportAngles
    vol.cellData.Theta = theta;
    vol.cellData.Phi = phi;
end
if cfg.exportDebugAngle
    % compute the absolute angle between the fiber direction and
    % an effecive long axis direction
    [~,~,V] = svd(gradAb, 'econ');
    longAxis = V(:,1)';
    vol.cellData.debugAngle = 90-rad2deg(acos(abs(fiber*longAxis')));
end

if cfg.exportFinalResult
    fprintf('Exporting result...                         '); tic;
    vtkWrite(vol, sprintf('%s.vtu', cfg.targetPrefix));
    fprintf('%.1f seconds\n', toc);
end

end