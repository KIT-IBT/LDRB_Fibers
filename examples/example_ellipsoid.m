addpath('../dependencies/gptoolbox/matrix');
addpath('../dependencies/gptoolbox/mesh');
addpath('../dependencies/gptoolbox/quat');
addpath('../dependencies/vtkToolbox/MATLAB');
addpath('../functions');

%%
alphaEndo = 60;
alphaEpi  = -60;
betaEndo  = 0;
betaEpi   = 0;

if ~exist('result_ellipsoid', 'dir')
    mkdir('result_ellipsoid');
end

clear cfg;
cfg.sourceDir = 'input_ellipsoid';
cfg.targetPrefix = 'result_ellipsoid/heart';
cfg.onlyOneVentricle = true;
cfg.volName = 'heart';
cfg.surNames = {'epi','lv','base','apex'};
cfg.alphaSeptLeft  = alphaEndo;
cfg.alphaSeptRight = alphaEndo;
cfg.alphaWallEndo  = alphaEndo;
cfg.alphaWallEpi   = alphaEpi;
cfg.betaSeptLeft   = betaEndo;
cfg.betaSeptRight  = betaEndo;
cfg.betaWallEndo   = betaEndo;
cfg.betaWallEpi    = betaEpi;
cfg.exportIntermediateResults = false;
cfg.exportFinalResult = true;
cfg.exportFiber = true;
cfg.exportSheet = true;
cfg.exportSheetnormal = true;
cfg.exportAngles = false;
cfg.outputAngleUnit = 'rad';
cfg.exportDebugAngle = true;
cfg.tol = 1e-12;
cfg.maxit = 1000;

% res = ldrb_main_original(cfg);
res = ldrb_main_adapted(cfg);
