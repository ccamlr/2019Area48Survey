%%
% Convert EK80 raw files with complex sample data into raw files with
% amplitude/angle data (a la EK60) and convert LSSS work files into
% Echoview region and line files. 

do_define_directories
baseDir = fullfile(baseDir, 'cruise_folders');

clear opt
% Cabo de Hornos
i = 1;
opt(i).workDir   = fullfile(baseDir, 'CDH2019828\ACOUSTIC\WORK');
opt(i).evDir     = fullfile(baseDir, 'CDH2019828\ACOUSTIC\EVL_EVR');
opt(i).evOutDir  = fullfile(baseDir, 'CDH2019828\ACOUSTIC');
opt(i).rawDir    = fullfile(baseDir, 'CDH2019828\EK80_RAW'); % WRONG dir...
opt(i).newRawDir = fullfile(baseDir, 'CDH2019828\ACOUSTIC\EK80\EK80_RAW_REDUCED');
opt(i).ship = 'CDH';
i = i + 1;

% Ukranian ship
opt(i).workDir   = fullfile(baseDir, 'MS2019999\ACOUSTIC\WORK');
opt(i).evDir     = fullfile(baseDir, 'MS2019999\ACOUSTIC\EVL_EVR');
opt(i).evOutDir  = fullfile(baseDir, 'MS2019999\ACOUSTIC');
opt(i).rawDir    = fullfile(baseDir, 'MS2019999\ACOUSTIC\ES80\ES80_ORIGINALRAWDATA');
opt(i).newRawDir = fullfile(baseDir, 'MS2019999\ACOUSTIC\ES80\ES80_RAWDATA_REDUCED');
opt(i).ship = 'MS';
i = i + 1;

% Kronprins Haakon
opt(i).workDir   = fullfile(baseDir, 'S2019701_PKRONPRINSHAAKON_9566\ACOUSTIC\LSSS\WORK');
opt(i).evDir     = fullfile(baseDir, 'S2019701_PKRONPRINSHAAKON_9566\ACOUSTIC\LSSS\EVL_EVR');
opt(i).evOutDir  = fullfile(baseDir, 'S2019701_PKRONPRINSHAAKON_9566\ACOUSTIC\LSSS');
opt(i).rawDir    = []; %fullfile(baseDir, 'S2019701_PKRONPRINSHAAKON_9566\ACOUSTIC\EK80\EK80_ORIGINALRAWDATA\EK80DK');
opt(i).newRawDir = fullfile(baseDir, 'S2019701_PKRONPRINSHAAKON_9566\ACOUSTIC\EK80\EK80_ORIGINALRAWDATA\EK80DK');
opt(i).ship = 'KPH';
i = i + 1;

% Kwang Ja Ho. EK60 data, so don't need to do EK80 data reduction.
opt(i).workDir   = fullfile(baseDir, 'KJH2019999\ACOUSTIC\LSSS_FILES\WORK');
opt(i).evDir     = fullfile(baseDir, 'KJH2019999\ACOUSTIC\LSSS_FILES\ECHOVIEW');
opt(i).evOutDir  = fullfile(baseDir, 'KJH2019999\ACOUSTIC\LSSS_FILES');
opt(i).rawDir    = [];
opt(i).newRawDir = fullfile(baseDir, 'KJH2019999\ACOUSTIC\EK60\EK60_ORIGINALRAWDATA');
opt(i).ship = 'KJH';
i = i + 1;

for i = 1:length(opt)
    opt(i).echoIntDir = fullfile(baseDir, 'data', 'echo-integration', opt(i).ship);
end


%% Convert EK80 full resolution data into smaller files

for j = 1:length(opt)

    if ~isempty(opt(j).rawDir) % only where there a rawDir
        disp(['Processing ' opt(j).ship])
        d = dir(fullfile(opt(j).rawDir, '*.raw'));
        
        skip = true;
        
        for i = 1:length(d)
            
            if skip && exist(fullfile(opt(j).newRawDir, d(i).name), 'file')
                disp(['Skipping ' d(i).name ' (already exists)'])
            else
                disp(['Squashing ' d(i).name ' (' num2str(i) ' of ' num2str(length(d)) ')'])
                EK80modifyFiles(fullfile(d(i).folder, d(i).name), ...
                    opt(j).newRawDir, 'sumTransducerChannels', true)
            end
        end
    end
end

%%
% Convert LSSS work files into Echoview .evr and .evl files
for k = 1:length(opt)
    disp(['Processing ' opt(k).ship])
    d = dir(fullfile(opt(k).newRawDir, '*.raw'));
    
    for i = 1:length(d) % for each .work file
        wFile = [d(i).name(1:end-3) 'work'];
        if exist(fullfile(opt(k).workDir, wFile), 'file')
            disp(['Converting ' wFile ' (' num2str(i) ' of ' num2str(length(d)) ')'])
            r = convertWorkToEchoview(wFile, opt(k).workDir, ...
                opt(k).newRawDir, opt(k).evDir, '120', ...
                'exportSchools', false, 'exportLayers', false);
            if r == 0
                disp(' No regions written')
            end
        else
            disp(['No work file found for ' d(i).name])
        end
    end
    
    % and merge these into 1 evl and 1 evr file for ease of importing into
    % Echoview

    % EVL files
    d = dir(fullfile(opt(k).evDir, '*.evl'));
    [~, ind] = sort({d.name});
    d = d(ind);
    
    j = 1;
    clear l
    for i = 1:length(d)
        %disp([num2str(i) ' of ' num2str(length(d))])
        fid = fopen(fullfile(d(i).folder, d(i).name), 'r');
        fgetl(fid);
        num = fscanf(fid, '%d\r\n', 1);
        for i = 1:num
            l{j} = fgetl(fid);
            j = j + 1;
        end
        fclose(fid);
    end
    
    evlFilename = fullfile(opt(k).evOutDir, [opt(k).ship '.evl']);
    fid = fopen(evlFilename, 'w');
    fprintf(fid, '%s\r\n', 'EVBD 3 3.00.41');
    fprintf(fid, '%d\r\n', length(l));
    for i = 1:length(l)
        fprintf(fid, '%s\r\n', l{i});
    end
    fclose(fid);
    disp(['Merged ' num2str(length(d)) ' .evl files'])
    
    % EVR files
    d = dir(fullfile(opt(k).evDir, '*.evr'));
    [~, ind] = sort({d.name});
    d = d(ind);
    
    j = 1;
    num = 0;
    clear l
    for i = 1:length(d)
        %disp([num2str(i) ' of ' num2str(length(d))])
        fid = fopen(fullfile(d(i).folder, d(i).name), 'r');
        fgetl(fid);
        n = fscanf(fid, '%d', 1);
        fgetl(fid);
        num = num + n;
        while ~feof(fid)
            l{j} = fgetl(fid);
            j = j + 1;
        end
        fclose(fid);
    end
    
    evrFilename = fullfile(opt(k).evOutDir, [opt(k).ship '.evr']);
    fid = fopen(evrFilename, 'w');
    fprintf(fid, '%s\r\n', 'EVRG 6 3.00.41');
    fprintf(fid, '%d\r\n', num);
    for i = 1:length(l)
        fprintf(fid, '%s\r\n', l{i});
    end
    fclose(fid);
    disp(['Merged ' num2str(length(d)) ' .evr files'])
    
    % and move the .evl and .evr files to the echo-integration directory 
    [~, n, e] = fileparts(evrFilename);
    movefile(evrFilename, fullfile(opt(k).echoIntDir, [n e]))
    
    [~, n, e] = fileparts(evlFilename);
    movefile(evlFilename, fullfile(opt(k).echoIntDir, [n e]))
end