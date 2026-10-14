% Generate moving_path.txt for orbit + translation + spin motion.
% File format used by the C++ code:
% time_s x y z qw qx qy qz
%
% Convention:
% - Position x/y/z is the rigid-body reference origin in global coordinates.
% - Quaternion qw/qx/qy/qz is body -> global.
% - The final pose combines:
%   1) orbit of the target reference origin around orbitCenterGlobal,
%   2) uniform translation in the global frame,
%   3) spin attitude about either a body-fixed axis or a global axis.

clear; clc;

outFile = 'moving_path_orbit_spin_translate.txt';

% ---- Time parameters ----
dt = 0.05;                 % s
timeLength = 36.0;         % s

% ---- Orbit parameters ----
orbitAxisName = 'z';       % global orbit axis: 'x', 'y', or 'z'
orbitOmega = 10.0;         % rad/s
orbitCenterGlobal = [0 0 0];

% Initial target reference-origin location relative to orbitCenterGlobal.
% Example: a target initially located on the negative Y axis, 2 m away.
targetOriginRel0 = [0 -2 0];

% ---- Translation parameters ----
velocityGlobal = [0 0 1];  % m/s, global translation velocity

% ---- Spin parameters ----
spinAxisName = 'x';        % spin axis: 'x', 'y', or 'z'
spinOmega = 0.0;           % rad/s. Set nonzero to enable spin.
spinAxisFrame = 'body';    % 'body' for body-fixed spin, 'global' for global-axis spin

% ---- Optional trajectory check points in the target body frame ----
% These are not written into moving_path.txt. They are only used to check
% the global trajectory of representative points on the target.
bodyPoints = [
    0 0 0
    1 0 0
    0 1 0
    0 0 1
];

% ---- Generate pose samples ----
t = (0:dt:timeLength).';
sampleCount = numel(t);

orbitAxisUnit = axisNameToUnit(orbitAxisName);
spinAxisUnit = axisNameToUnit(spinAxisName);

orbitTheta = orbitOmega .* t;
spinTheta = spinOmega .* t;

qOrbit = axisAngleToQuatRows(orbitAxisUnit, orbitTheta);
qSpin = axisAngleToQuatRows(spinAxisUnit, spinTheta);

pos = zeros(sampleCount, 3);
q = zeros(sampleCount, 4);

for i = 1:sampleCount
    Rorbit = quatToRotmBodyToGlobal(qOrbit(i, :));
    pos(i, :) = orbitCenterGlobal + targetOriginRel0 * Rorbit.' + t(i) .* velocityGlobal;

    switch lower(spinAxisFrame)
        case 'body'
            % Body-fixed spin: R_total = R_orbit * R_spin.
            q(i, :) = quatMultiply(qOrbit(i, :), qSpin(i, :));
        case 'global'
            % Global-axis spin: R_total = R_spin * R_orbit.
            q(i, :) = quatMultiply(qSpin(i, :), qOrbit(i, :));
        otherwise
            error('spinAxisFrame must be body or global.');
    end
    q(i, :) = q(i, :) ./ norm(q(i, :));
end

q(abs(q) < 1e-15) = 0;
pos(abs(pos) < 1e-15) = 0;

data = [t pos q];
data(abs(data) < 1e-15) = 0;

fid = fopen(outFile, 'w');
if fid < 0
    error('Cannot open %s for writing.', outFile);
end

fprintf(fid, '# time_s x y z qw qx qy qz\n');
fprintf(fid, '# Orbit: center = (%.12g, %.12g, %.12g), axis = +%s, omega = %.12g rad/s.\n', ...
    orbitCenterGlobal(1), orbitCenterGlobal(2), orbitCenterGlobal(3), ...
    upper(orbitAxisName), orbitOmega);
fprintf(fid, '# Initial target origin relative to orbit center = (%.12g, %.12g, %.12g).\n', ...
    targetOriginRel0(1), targetOriginRel0(2), targetOriginRel0(3));
fprintf(fid, '# Translation velocity = (%.12g, %.12g, %.12g) m/s.\n', ...
    velocityGlobal(1), velocityGlobal(2), velocityGlobal(3));
fprintf(fid, '# Spin: axis = +%s, omega = %.12g rad/s, frame = %s.\n', ...
    upper(spinAxisName), spinOmega, spinAxisFrame);
fprintf(fid, '# Quaternion is body -> global. dt = %.12g s, duration = %.12g s.\n', ...
    dt, timeLength);
fprintf(fid, '%.6f %.12f %.12f %.12f %.12f %.12f %.12f %.12f\n', data.');
fclose(fid);

fprintf('Wrote %d states to %s\n', size(data, 1), outFile);

% ---- Optional check: global coordinates of body reference points ----
pointTraj = zeros(sampleCount, size(bodyPoints, 1), 3);
for i = 1:sampleCount
    R = quatToRotmBodyToGlobal(q(i, :));
    pointTraj(i, :, :) = bodyPoints * R.' + pos(i, :);
end

disp('First-time body reference point positions:');
disp(squeeze(pointTraj(1, :, :)));
disp('Last-time body reference point positions:');
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

function q = axisAngleToQuatRows(axisUnit, theta)
q = zeros(numel(theta), 4);
q(:, 1) = cos(theta ./ 2);
q(:, 2:4) = sin(theta ./ 2) .* axisUnit;
end

function q = quatMultiply(a, b)
aw = a(1); ax = a(2); ay = a(3); az = a(4);
bw = b(1); bx = b(2); by = b(3); bz = b(4);

q = [
    aw*bw - ax*bx - ay*by - az*bz, ...
    aw*bx + ax*bw + ay*bz - az*by, ...
    aw*by - ax*bz + ay*bw + az*bx, ...
    aw*bz + ax*by - ay*bx + az*bw
];
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
