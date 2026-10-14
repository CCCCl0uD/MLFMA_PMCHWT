% Generate moving_path.txt for FSB/HSB-MLFMA-Moving and ACA-RHS-MLFMA-Moving.
% File format used by the C++ code:
% time_s x y z qw qx qy qz
%
% The quaternion is body -> global. The moving algorithms then use the
% equivalent monostatic direction transform, instead of transforming the
% target impedance matrix.

clear; clc;

outFile = 'moving_path_QiuZhui_z.txt';

% ---- User parameters ----
axisName = 'z';              % 'x', 'y', or 'z'
omega = 10.0;                % rad/s
dt = 0.05;                   % s
timeLength = 36.0;           % s
targetOriginGlobal0 = [0 0 0]; % initial rigid-body origin/centroid position
velocityGlobal = [0 0 0];      % m/s, translational velocity in global frame

% Reference locations of the four cone-spheres in the body frame.
% These are not written into moving_path.txt; they are here for optional
% trajectory checking because the mesh already contains this geometry.
bodyPoints = [
    0 -2 0
];

% ---- Generate pose samples ----
t = (0:dt:timeLength).';
theta = omega .* t;

axisUnit = axisNameToUnit(axisName);
q = zeros(numel(t), 4);
q(:, 1) = cos(theta ./ 2);
q(:, 2:4) = sin(theta ./ 2) .* axisUnit;
q(abs(q) < 1e-15) = 0;

pos = targetOriginGlobal0 + t .* velocityGlobal;
data = [t pos q];
data(abs(data) < 1e-15) = 0;

fid = fopen(outFile, 'w');
if fid < 0
    error('Cannot open %s for writing.', outFile);
end

fprintf(fid, '# time_s x y z qw qx qy qz\n');
fprintf(fid, '# Rigid body rotation about +%s axis, omega = %.12g rad/s, dt = %.12g s, duration = %.12g s.\n', ...
    upper(axisName), omega, dt, timeLength);
fprintf(fid, '# Quaternion is body -> global. Target origin starts at (%.12g, %.12g, %.12g), velocity = (%.12g, %.12g, %.12g) m/s.\n', ...
    targetOriginGlobal0(1), targetOriginGlobal0(2), targetOriginGlobal0(3), ...
    velocityGlobal(1), velocityGlobal(2), velocityGlobal(3));
fprintf(fid, '%.6f %.12f %.12f %.12f %.12f %.12f %.12f %.12f\n', data.');
fclose(fid);

fprintf('Wrote %d states to %s\n', size(data, 1), outFile);

% Optional check: global coordinates of the four body reference points.
% pointTraj(i, j, :) is point j at time index i in global coordinates.
pointTraj = zeros(numel(t), size(bodyPoints, 1), 3);
for i = 1:numel(t)
    R = quatToRotmBodyToGlobal(q(i, :));
    pointTraj(i, :, :) = bodyPoints * R.' + pos(i, :);
end

disp('First-time reference point positions:');
disp(squeeze(pointTraj(1, :, :)));
disp('Last-time reference point positions:');
disp(squeeze(pointTraj(end, :, :)));

function axisUnit = axisNameToUnit(axisName)
switch lower(axisName)
    case 'x'
        axisUnit = [1 0 0];
    case 'y'
        axisUnit = [0 1 0];
    case 'z'
        axisUnit = [0 0 1];
    otherwise
        error('axisName must be x, y, or z.');
end
end

function R = quatToRotmBodyToGlobal(q)
q = q ./ norm(q);
w = q(1); x = q(2); y = q(3); z = q(4);
R = [
    1 - 2 * (y * y + z * z), 2 * (x * y - z * w),     2 * (x * z + y * w)
    2 * (x * y + z * w),     1 - 2 * (x * x + z * z), 2 * (y * z - x * w)
    2 * (x * z - y * w),     2 * (y * z + x * w),     1 - 2 * (x * x + y * y)
];
end
