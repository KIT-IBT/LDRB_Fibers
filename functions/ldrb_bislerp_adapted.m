function Qab = ldrb_bislerp_adapted(Qa, Qb, t)
% Function 4 in supplement of Bayer 2012, but without quaternion 'flipping':
% https://doi.org/10.1007/s10439-012-0593-5

qa = rotm2quat(Qa);
qb = rotm2quat(Qb);

t = max(min(t,1),0);

% qab = quatinterp(qa, qb, t, 'slerp'); % requires MATLAB's Aerospace Toolbox
qab = quatslerp(qa, qb, t);             % uses gptoolbox instead

Qab = quat2rotm(qab);  % requires MATLAB's Robotics System Toolbox
% Qab = quat2mat(qab); % uses gptoolbox instead

end