% Merge the last data line of FEKO .ffe files into one TXT file.
% Files are sorted by the number in parentheses, for example:
%   JL_1e9_FarField1(1).ffe, JL_1e9_FarField1(2).ffe, ..., JL_1e9_FarField1(10).ffe

clear; clc;

inputDir = 'D:\MyCode\PMCHWT_MLFMA\DATA\JL_1e9\feko';
filePrefix = 'JL_1e9_FarField1';
outputFile = fullfile(inputDir, 'JL_1e9_FarField1_last_lines.txt');

files = dir(fullfile(inputDir, [filePrefix, '(*).ffe']));
if isempty(files)
    error('No .ffe files found in: %s', inputDir);
end

fileNumbers = nan(numel(files), 1);
for i = 1:numel(files)
    token = regexp(files(i).name, [regexptranslate('escape', filePrefix), '\((\d+)\)\.ffe$'], ...
                   'tokens', 'once');
    if isempty(token)
        error('Unexpected file name: %s', files(i).name);
    end
    fileNumbers(i) = str2double(token{1});
end

[fileNumbers, sortIdx] = sort(fileNumbers);
files = files(sortIdx);

outFid = fopen(outputFile, 'w');
if outFid < 0
    error('Cannot open output file: %s', outputFile);
end
cleanupObj = onCleanup(@() fclose(outFid));

fprintf(outFid, ['Index\tTheta\tPhi\tRe_Etheta\tIm_Etheta\tRe_Ephi\tIm_Ephi\t', ...
                 'RCS_Theta\tRCS_Phi\tRCS_Total\n']);

for i = 1:numel(files)
    filePath = fullfile(files(i).folder, files(i).name);
    lastLine = readLastNonemptyLine(filePath);
    values = sscanf(lastLine, '%f').';

    if numel(values) ~= 9
        error('Expected 9 numeric values in the last line of %s, but found %d.', ...
              files(i).name, numel(values));
    end
    
    fprintf(outFid, '%d', fileNumbers(i));
    fprintf(outFid, '\t%.8E', values);
    fprintf(outFid, '\n');

end

fprintf('Merged %d files into:\n%s\n', numel(files), outputFile);

function lastLine = readLastNonemptyLine(filePath)
    fid = fopen(filePath, 'r');
    if fid < 0
        error('Cannot open file: %s', filePath);
    end
    cleanupObj = onCleanup(@() fclose(fid));

    lastLine = '';
    while true
        line = fgetl(fid);
        if ~ischar(line)
            break;
        end
        if ~isempty(strtrim(line))
            lastLine = line;
        end
    end
    
    if isempty(lastLine)
        error('File is empty: %s', filePath);
    end

end