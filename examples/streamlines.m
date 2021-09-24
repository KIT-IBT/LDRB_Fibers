addpath('../dependencies/vtkToolbox/MATLAB');

%% Visualize fibers by computing streamlines

N = 5000; % number of streamlines

v = vtkRead('result_coarse/heart.vtu');
v = vtkDeleteDataArrays(v, {'cellData.Fiber'});
seeds = vtkCellCentroids(v);
rng(1);
seedIds = randi(size(seeds.points,1), N, 1);
seeds = vtkCreateStruct(seeds.points(seedIds,:));
edgLen = mean(vtkEdgeLengths(v));
sl = vtkStreamTracer(v, 'cells', 'Fiber', seeds, 'both', edgLen);
vtkWrite(sl, 'result_coarse/streamlines.vtp'); 

% Open streamlines.vtp in ParaView and apply a 'Tube' filter
