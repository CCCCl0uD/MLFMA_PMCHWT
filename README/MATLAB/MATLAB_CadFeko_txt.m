clear; clc;

cadfeko = '"C:\Program Files\Altair\2022\feko\bin\cadfeko.exe"';

workDir = 'D:\MyCode\PMCHWT_MLFMA\DATA\Cube_0d5m_Die_1e9\feko';
baseName = 'Cube_0d5m_Die_1e9';

cfxFile = fullfile(workDir, [baseName '.cfx']);
templateLua = fullfile(workDir, 'run_cadfeko.lua');
outFile = fullfile(workDir, [baseName '.out']);

cases = {
    'eps2d25_tand0d01', '0'
    'eps2d25_tand0d01', '90'
    'eps4d0_tand0d02',  '0'
    'eps4d0_tand0d02',  '90'
    'eps10d0_tand0d05', '0'
    'eps10d0_tand0d05', '90'
};

for i = 1:size(cases,1)
    mediumName = cases{i,1};
    polAngle = cases{i,2};

    if strcmp(polAngle, '0')
        polType = 'VV';
        fieldComponent = 'Etheta';
    elseif strcmp(polAngle, '90')
        polType = 'HH';
        fieldComponent = 'Ephi';
    else
        polType = ['POL' polAngle];
        fieldComponent = 'Ephi';
    end

    luaFile = fullfile(workDir, sprintf('run_%s_pol%s.lua', mediumName, polAngle));

    datFile = fullfile(workDir, sprintf('%s_%s_%s.dat', ...
        baseName, mediumName, polType));

    txtFile = fullfile(workDir, sprintf('%s_%s_%s_%s.txt', ...
        baseName, mediumName, polAngle, polType));

    txt = fileread(templateLua);
    txt = strrep(txt, '**MEDIUM_NAME**', mediumName);
    txt = strrep(txt, '**POL_ANGLE**', polAngle);

    fid = fopen(luaFile, 'w');
    if fid < 0
        error('Cannot open Lua file for writing: %s', luaFile);
    end
    fwrite(fid, txt);
    fclose(fid);

    if isfile(outFile)
        delete(outFile);
    end

    fprintf('\nRunning: %s pol %s (%s)\n', mediumName, polAngle, polType);

    cmd = sprintf('%s --run-script "%s" "%s"', cadfeko, luaFile, cfxFile);
    [status, msg] = system(cmd);

    if ~isfile(outFile)
        fprintf('Failed: %s pol %s\n', mediumName, polAngle);
        fprintf('CADFEKO exit status: %d\n', status);
        fprintf('CADFEKO message:\n%s\n', msg);
        continue;
    end

    try
        nPoints = writePolarizedRcsFromOut(outFile, datFile, txtFile, fieldComponent);

        fprintf('Finished: %s pol %s -> %s (%d points)\n', ...
            mediumName, polAngle, datFile, nPoints);

        fprintf('TXT saved -> %s\n', txtFile);
    catch ME
        fprintf('RCS extraction failed: %s pol %s\n', mediumName, polAngle);
        fprintf('%s\n', ME.message);
    end
end

function nPoints = writePolarizedRcsFromOut(outFile, datFile, txtFile, fieldComponent)
    lines = splitlines(fileread(outFile));

    thetaList = [];
    phiList = [];
    rcsDbsmList = [];

    inFarFieldTable = false;

    for k = 1:numel(lines)
        line = strtrim(lines{k});

        if contains(line, 'LOCATION') && contains(line, 'ETHETA') && contains(line, 'EPHI')
            inFarFieldTable = true;
            continue;
        end

        if ~inFarFieldTable
            continue;
        end

        nums = sscanf(line, '%f');

        % FEKO far-field columns:
        % 1 theta
        % 2 phi
        % 3 |Etheta|
        % 4 phase(Etheta)
        % 5 |Ephi|
        % 6 phase(Ephi)
        % 7 total scattering cross section
        if numel(nums) >= 7
            theta = nums(1);
            phi = nums(2);

            if strcmpi(fieldComponent, 'Etheta')
                esMag = nums(3);
            elseif strcmpi(fieldComponent, 'Ephi')
                esMag = nums(5);
            else
                error('Unknown field component: %s', fieldComponent);
            end

            rcsDbsm = 10 * log10(4*pi*abs(esMag)^2);

            thetaList(end+1,1) = theta;
            phiList(end+1,1) = phi;
            rcsDbsmList(end+1,1) = rcsDbsm;
        elseif ~isempty(phiList)
            break;
        end
    end

    nPoints = numel(phiList);

    if nPoints == 0
        error('No far-field data rows were found in: %s', outFile);
    end

    fid = fopen(datFile, 'w');
    if fid < 0
        error('Cannot open DAT file for writing: %s', datFile);
    end

    fprintf(fid, 'Far field\n');
    fprintf(fid, '"Phi"[deg]\t"FarField1 [dBsm]"\t\n');

    for k = 1:nPoints
        fprintf(fid, '%.8E\t%.8E\t\n', phiList(k), rcsDbsmList(k));
    end

    fclose(fid);

    fid = fopen(txtFile, 'w');
    if fid < 0
        error('Cannot open TXT file for writing: %s', txtFile);
    end

    fprintf(fid, 'Theta[deg]\tPhi[deg]\tRCS[dBsm]\n');

    for k = 1:nPoints
        fprintf(fid, '%.8E\t%.8E\t%.8E\n', ...
            thetaList(k), phiList(k), rcsDbsmList(k));
    end

    fclose(fid);
end