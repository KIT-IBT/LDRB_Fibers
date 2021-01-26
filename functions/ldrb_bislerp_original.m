function Qab = ldrb_bislerp_original(Qa, Qb, t)
% Function 4 in supplement of Bayer 2012:
% https://doi.org/10.1007/s10439-012-0593-5

qa = rotm2quat(Qa);
qb = rotm2quat(Qb);

qi = [ ...
    1 0 0 0; -1 0 0 0; ...
    0 1 0 0; 0 -1 0 0; ...
    0 0 1 0; 0 0 -1 0; ...
    0 0 0 1; 0 0 0 -1  ...
    ];
qm = quatmultiply(repmat(qa,8,1), qi);
dot_qm_qb = dot(qm, repmat(qb,8,1), 2);
[~,ind] = max(abs(dot_qm_qb));
qm = qm(ind,:);

t = max(min(t,1),0);

% qab = quatinterp(qm, qb, t, 'slerp'); % requires MATLAB's Aerospace Toolbox
qab = quatslerp(qm, qb, t);             % uses gptoolbox instead

Qab = quat2rotm(qab);  % requires MATLAB's Robotics System Toolbox
% Qab = quat2mat(qab); % uses gptoolbox instead

end