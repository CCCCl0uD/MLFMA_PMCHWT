
endclear; clc;

cadfeko = '"C:\Program Files\Altair\2022\feko\bin\cadfeko.exe"';

workDir = 'D:\MyCode\PMCHWT_MLFMA\DATA\Cube_0d5m_Die_1e9\feko';
baseName = 'Cube_0d5m_Die_1e9';

cfxFile = fullfile(workDir, [baseName '.cfx']);
templateLua = fullfile(workDir, 'run_cadfeko.lua');

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
    elseif strcmp(polAngle, '90')
        polType = 'HH';
    else
        polType = ['POL' polAngle];
    end

    luaFile = fullfile(workDir, sprintf('run_%s_pol%s.lua', mediumName, polAngle));

    txt = fileread(templateLua);
    txt = strrep(txt, '**MEDIUM_NAME**', mediumName);
    txt = strrep(txt, '**POL_ANGLE**', polAngle);

    fid = fopen(luaFile, 'w');
    fwrite(fid, txt);
    fclose(fid);

    srcFek = fullfile(workDir, [baseName '.fek']);
    dstFek = fullfile(workDir, sprintf('%s_%s_%s.fek', baseName, mediumName, polType));

    if isfile(srcFek)
        delete(srcFek);
    end

    fprintf('\nRunning: %s pol %s\n', mediumName, polAngle);

    cmd = sprintf('%s --run-script "%s" "%s"', cadfeko, luaFile, cfxFile);
    [status, msg] = system(cmd);

    if isfile(srcFek)
        copyfile(srcFek, dstFek, 'f');
        fprintf("Finished: %s pol %s -> %s\n", mediumName, polAngle, dstFek);
    else
        fprintf("Failed: %s pol %s\n", mediumName, polAngle);
        fprintf("CADFEKO exit status: %d\n", status);
        fprintf("CADFEKO message:\n%s\n", msg);
    end