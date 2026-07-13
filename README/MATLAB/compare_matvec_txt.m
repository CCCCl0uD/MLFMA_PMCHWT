%% compare_matvec_txt.m
% Compare two one-shot matrix-vector multiplication txt files.
% Expected columns:
% index x_real x_imag y_real y_imag abs_y

clear;
clc;

%% User options
showFigure = true;
topK = 10;

%% Select two txt files
[fileA, pathA] = uigetfile({'*.txt', 'Text files (*.txt)'; '*.*', 'All files'}, ...
    'Select first matvec txt file');
if isequal(fileA, 0)
    error('No first file selected.');
end

[fileB, pathB] = uigetfile({'*.txt', 'Text files (*.txt)'; '*.*', 'All files'}, ...
    'Select second matvec txt file');
if isequal(fileB, 0)
    error('No second file selected.');
end

fullA = fullfile(pathA, fileA);
fullB = fullfile(pathB, fileB);

fprintf('File A: %s\n', fullA);
fprintf('File B: %s\n\n', fullB);

%% Read data
dataA = readMatvecTxt(fullA);
dataB = readMatvecTxt(fullB);

%% Basic checks
nA = height(dataA);
nB = height(dataB);

fprintf('Rows in A: %d\n', nA);
fprintf('Rows in B: %d\n', nB);

if nA ~= nB
    error('Row count mismatch: A has %d rows, B has %d rows.', nA, nB);
end

if any(dataA.index ~= dataB.index)
    warning('Index columns are not identical. Comparison will still use row order.');
end

%% Extract columns
yA = complex(dataA.y_real, dataA.y_imag);
yB = complex(dataB.y_real, dataB.y_imag);

absA = dataA.abs_y;
absB = dataB.abs_y;

diffY = yB - yA;
diffReal = dataB.y_real - dataA.y_real;
diffImag = dataB.y_imag - dataA.y_imag;
diffAbs = absB - absA;

absDiffY = abs(diffY);
absDiffReal = abs(diffReal);
absDiffImag = abs(diffImag);
absDiffAbs = abs(diffAbs);

%% Error metrics
epsDen = 1.0e-300;

maxAbsErrY = max(absDiffY);
relL2ErrY = norm(diffY, 2) / max(norm(yA, 2), epsDen);

maxAbsErrReal = max(absDiffReal);
relL2ErrReal = norm(diffReal, 2) / max(norm(dataA.y_real, 2), epsDen);

maxAbsErrImag = max(absDiffImag);
relL2ErrImag = norm(diffImag, 2) / max(norm(dataA.y_imag, 2), epsDen);

maxAbsErrAbs = max(absDiffAbs);
relL2ErrAbs = norm(diffAbs, 2) / max(norm(absA, 2), epsDen);

fprintf('\n===== Complex y comparison =====\n');
fprintf('max abs error of y        : %.17e\n', maxAbsErrY);
fprintf('relative L2 error of y    : %.17e\n', relL2ErrY);

fprintf('\n===== Column comparison =====\n');
fprintf('y_real max abs error      : %.17e\n', maxAbsErrReal);
fprintf('y_real relative L2 error  : %.17e\n', relL2ErrReal);
fprintf('y_imag max abs error      : %.17e\n', maxAbsErrImag);
fprintf('y_imag relative L2 error  : %.17e\n', relL2ErrImag);
fprintf('abs_y  max abs error      : %.17e\n', maxAbsErrAbs);
fprintf('abs_y  relative L2 error  : %.17e\n', relL2ErrAbs);

%% Show largest error rows
[sortedErr, order] = sort(absDiffY, 'descend');
topK = min(topK, numel(order));
topRows = order(1:topK);

fprintf('\n===== Top %d rows by abs(yB - yA) =====\n', topK);
fprintf('row\tindex\tyA_real\tyA_imag\tyB_real\tyB_imag\tabs_err_y\tabs_err_abs_y\n');

for k = 1:topK
    r = topRows(k);
    fprintf('%d\t%d\t%.9e\t%.9e\t%.9e\t%.9e\t%.9e\t%.9e\n', ...
        r, dataA.index(r), ...
        dataA.y_real(r), dataA.y_imag(r), ...
        dataB.y_real(r), dataB.y_imag(r), ...
        sortedErr(k), absDiffAbs(r));
end

%% Optional figures
if showFigure
    figure('Name', 'MatVec txt comparison', 'Color', 'w');

    subplot(3, 1, 1);
    semilogy(dataA.index, absDiffY, 'LineWidth', 1.2);
    grid on;
    xlabel('index');
    ylabel('|y_B - y_A|');
    title('Complex y absolute error');

    subplot(3, 1, 2);
    semilogy(dataA.index, absDiffReal, 'LineWidth', 1.2);
    hold on;
    semilogy(dataA.index, absDiffImag, 'LineWidth', 1.2);
    grid on;
    xlabel('index');
    ylabel('absolute error');
    legend('|real diff|', '|imag diff|', 'Location', 'best');
    title('Real and imaginary column errors');

    subplot(3, 1, 3);
    semilogy(dataA.index, absDiffAbs, 'LineWidth', 1.2);
    grid on;
    xlabel('index');
    ylabel('|abs\_y_B - abs\_y_A|');
    title('abs\_y column error');
end

%% Local function
function data = readMatvecTxt(filePath)
    opts = detectImportOptions(filePath, ...
        'FileType', 'text', ...
        'CommentStyle', '#', ...
        'Delimiter', '\t');

    raw = readmatrix(filePath, ...
        'FileType', 'text', ...
        'CommentStyle', '#');

    if isempty(raw)
        error('File has no numeric data: %s', filePath);
    end

    if size(raw, 2) < 6
        error('Expected at least 6 numeric columns in file: %s', filePath);
    end

    data = table();
    data.index = raw(:, 1);
    data.x_real = raw(:, 2);
    data.x_imag = raw(:, 3);
    data.y_real = raw(:, 4);
    data.y_imag = raw(:, 5);
    data.abs_y = raw(:, 6);

    if any(~isfinite(data.y_real)) || any(~isfinite(data.y_imag)) || any(~isfinite(data.abs_y))
        error('File contains NaN or Inf in y_real, y_imag, or abs_y: %s', filePath);
    end
end