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

if ~exist('result_coarse', 'dir')
    mkdir('result_coarse');
end

clear cfg;
cfg.sourceDir = 'input_coarse';
cfg.targetPrefix = 'result_coarse/heart';
cfg.onlyOneVentricle = false;
cfg.volName = 'heart';
cfg.surNames = {'heart'};
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
