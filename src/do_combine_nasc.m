%%
% A script to combine the Echoview integrals from each vessel into one
% dataset, with each NASC value labelled with vessel and transect name.

do_define_directories
dataDir = fullfile(baseDir, 'data', 'echo-integration');


vessels = {'CDH' 'KJH' 'KPH' 'MS' 'FRH' 'RRS'};

nasc = table([]);
for i = 1:length(vessels)
    disp(['Loading ' vessels{i}])
    
    d = dir(fullfile(dataDir, vessels{i}, 'krillNASC*.csv'));
    
    for j = 1:length(d)
        if strcmp(vessels{i}, 'FRH')
            % format is different to the other vessels
            opts = delimitedTextImportOptions("NumVariables", 11);
            opts.DataLines = [2, Inf];
            opts.Delimiter = ",";
            opts.VariableNames = ["Interval", "Sv_mean", "NASC", "Date_S", "Time_S", "Date_E", "Time_E", "Date_M", "Time_M", "Lat_M", "Lon_M"];
            opts.VariableTypes = ["double", "double", "double", "datetime", "datetime", "datetime", "datetime", "datetime", "datetime", "double", "double"];
            opts = setvaropts(opts, 4, "InputFormat", "yyyyMMdd");
            opts = setvaropts(opts, 5, "InputFormat", "HH:mm:ss.SSSS");
            opts = setvaropts(opts, 6, "InputFormat", "yyyyMMdd");
            opts = setvaropts(opts, 7, "InputFormat", "HH:mm:ss.SSSS");
            opts = setvaropts(opts, 8, "InputFormat", "yyyyMMdd");
            opts = setvaropts(opts, 9, "InputFormat", "HH:mm:ss.SSSS");
            opts.ExtraColumnsRule = "ignore";
            opts.EmptyLineRule = "read";
            
            n = readtable(fullfile(d(j).folder, d(j).name), opts);
            % remove unwanted variables and rename some
            n = removevars(n, {'Interval', 'Sv_mean', 'Date_S', 'Time_S', 'Date_E', 'Time_E'});
            n.Properties.VariableNames{'Time_M'} = 'Ping_time';
            n.Properties.VariableNames{'Date_M'} = 'Ping_date';
            n.Properties.VariableNames{'Lat_M'} = 'Latitude';
            n.Properties.VariableNames{'Lon_M'} = 'Longitude';            

            % To get the transects, we divide on time, based on manual
            % inspection. We use the AMLR stratum names but make up
            % transect numbers (see EMM-11/26). Note that FRH did extra
            % transects that asked for. These extra ones are in the West
            % stratum and overlap with those done by KJH. We don't use
            % these extra transects.
            transects = {
                '20190205T102300', '20190205T153100', "Elephant", "02";
                '20190205T172300', '20190206T090049', "Elephant", "01";
                '20190206T092050', '20190206T164000', "Joinville", "01";
                '20190206T192400', '20190207T020500', "Joinville", "02";
                '20190207T052600', '20190207T101316', "Joinville", "03";
                '20190207T110955', '20190207T162900', "Bransfield", "01";
                '20190207T194500', '20190207T234100', "Bransfield", "02";
                '20190208T025800', '20190208T061900', "Bransfield", "03";
                '20190208T101200', '20190208T133300', "Bransfield", "04";
                '20190208T161900', '20190208T203000', "Bransfield", "05";
                '20190208T231100', '20190209T021800', "Bransfield", "06";
                '20190209T050200', '20190209T092301', "Bransfield", "07"};
                %'20190209T113440', '20190209T153450', "West", "07";
                %'20190209T180130', '20190209T221645', "West", "06";
                %'20190210T003000', '20190210T054150', "West", "05";
                %'20190210T081200', '20190210T135850', "West", "04";
                %'20190210T161030', '20190210T211230', "West", "03"};
                
            t = n.Ping_date + timeofday(n.Ping_time);
            for k = 1:size(transects,1)
                kk = find(t >= datetime(transects{k,1}) & t <= datetime(transects{k,2}));
                n.Transect(kk) = repmat(transects{k,4}, length(kk), 1);
                n.Stratum(kk) = transects{k,3};
            end
            % if data does not have a Transect or Stratum, then it is not
            % on a transect, so remove it
            k = ~ismissing(n.Transect);
            n = n(k,:);
        else % data done via the LSSS/Echoview path by IMR
            % Setup the parameters to import the Echoview output files. Echoview
            % doesn't provide a column label for the NASC column, so this mucks up
            % Matlab's automatic importing tools. So set all options explicitly and
            % give a name to the NASC column.
            opts = delimitedTextImportOptions("NumVariables", 14);
            opts.DataLines = [2, Inf];
            opts.Delimiter = ",";
            opts.VariableNames = ["Ping_index", "Distance_gps", "Distance_vl", "Ping_date", "Ping_time", "Ping_milliseconds", "Latitude", "Longitude", "Depth_start", "Depth_stop", "Range_start", "Range_stop", "Sample_count", "NASC"];
            opts.VariableTypes = ["double", "double", "double", "datetime", "datetime", "double", "double", "double", "double", "double", "double", "double", "double", "double"];
            opts = setvaropts(opts, 4, "InputFormat", "yyyy-MM-dd");
            opts = setvaropts(opts, 5, "InputFormat", "HH:mm:ss");
            opts.ExtraColumnsRule = "ignore";
            opts.EmptyLineRule = "read";
            
            n = readtable(fullfile(d(j).folder, d(j).name), opts);
            % reduce to minimum set of required variables
            n = removevars(n, {'Ping_index', 'Distance_gps', 'Distance_vl', 'Ping_milliseconds', 'Depth_start', 'Depth_stop', 'Range_start', 'Range_stop', 'Sample_count'});

            % Trim the NASC using the on/off transect timestamps given in
            % the echo-integration directories to ensure that we don't get
            % unwanted integration data. 
            % only do this if there is a 'transects.csv' file
            if exist(fullfile(d(j).folder, 'transects.csv'), 'file')
                transects = readtable(fullfile(d(j).folder, 'transects.csv'));
                % transects.csv files have the time with milliseconds and
                % timezone and readtable returns the time as a character
                % string, so deal with that (and truncates the milliseconds)
                transects.Time = datetime(transects.Time, 'InputFormat', 'HH:mm:ss.SSSZ', 'TimeZone', 'UTC');
                transects.Time.TimeZone = '';
                transects.Time.Second = floor(transects.Time.Second);
                
                % Get the timestamp that is in the echo integral filename
                parts = split(d(j).name, '_');
                onTime = datetime(parts{3});
                
                % Combine transect date and time into a timestamp
                transects.timestamp = transects.Date + timeofday(transects.Time);
                
                % find the on/off lines for the current echo integral .csv file
                tt = find(transects.timestamp == onTime & startsWith(transects.Event, "On"));
                if length(tt) ~= 1
                    tt
                    error('Error in do_combine_nasc')
                end
                onTime = transects.timestamp(tt);
                offTime = transects.timestamp(tt+1);
                
                % now trim the nasc based on the on and off time
                tt = n.Ping_date+timeofday(n.Ping_time);
                ii = find(tt >= onTime & tt <= offTime);
                n = n(ii,:);
            end
            
            % add transect and stratum id to the data
            t = split(d(j).name, '_');            
            [st, tr] = sortOutNaming(string(t{2}));
            
            n.Stratum = repmat(st, size(n.NASC, 1), 1);
            n.Transect = repmat(tr, size(n.NASC, 1), 1);
        end
        
        % add vessel name to the data        
        n.Vessel = repmat(string(vessels{i}), size(n.NASC));

        if isempty(nasc)
            nasc = n;
        else
            nasc = [nasc; n];
        end
        
    end
    
end

% Tidy up things:
%   For miles inside exclude regions Echoview produces values of -9.9e37.
%   Remove those rows.
nasc = nasc(nasc.NASC ~= -9.9e37,:);
% merge the date and time columns into one
nasc.Ping_timestamp = nasc.Ping_date + timeofday(nasc.Ping_time);
nasc.Ping_timestamp.Format = 'yyyy-MM-dd HH:mm:ss.SSSS';
nasc = removevars(nasc, {'Ping_date', 'Ping_time'});
nasc = movevars(nasc, {'Ping_timestamp'}, 'Before', 1);

% flag each nasc value to show if it was taken within (or not) civil
% daylight hours.
for i = 1:height(nasc)
    nasc.civilDaytime(i) = dayOrNight(nasc.Longitude(i), nasc.Latitude(i), nasc.Ping_timestamp(i));
end

% and save the combined dataset
save(fullfile(resultsDir, 'NASC - data'), 'nasc')


%%%%%%%%%%%%%%%%
% Maps

if false % these plots have been moved to showing rho, and are now down in the do_estimate_biomass.m script. 
    strata = jsondecode(fileread(fullfile(baseDir, repoDir, 'map_data', 'survey strata.geojson')));
    
    % Map coloured by vessel
    figure(1)
    clf
    plot_standard_map(strata)
    
    % Use a different colour for each vessel
    maxNASC = max(nasc.NASC);
    v = unique(nasc.Vessel);
    h = nan(length(v), 1);
    for i = 1:length(v)
        j = find(nasc.Vessel == v(i));
        h(i) = m_scatter(nasc.Longitude(j), nasc.Latitude(j), nasc.NASC(j)/maxNASC*200+1, 'filled');
    end
    
    m_grid('box', 'on')
    legend(h, v, 'Location', 'SouthEast')
    
    ifile = fullfile(resultsDir, 'NASC - by vessel.png');
    print(ifile, '-dpng','-r300')
    crop_image(ifile)
    
    % Map coloured by stratum. Do several, based on different areas
    % Entire area. Too crowded to be really useful...
    figure(2)
    clf
    plot_standard_map(strata)
    
    % Use a different colour and symbol for each statum
    maxNASC = max(nasc.NASC);
    s = unique(nasc.Stratum);
    symbols = {'o' 'o' 'o' 'o' 'o' 'o' 'o' 'd' 'd' 'd' 'd' 'd' 'd' 'd'};
    h = [];
    for i = 1:length(s)
        j = find(nasc.Stratum == s(i));
        h(i) = m_scatter(nasc.Longitude(j), nasc.Latitude(j), nasc.NASC(j)/maxNASC*200+1, 'filled', symbols{i});
    end
    
    m_grid('box', 'on')
    legend(h, s, 'Location', 'SouthEast', 'NumColumns', 2, 'Interpreter', 'none')
    
    ifile = fullfile(resultsDir, 'NASC - by stratum.png');
    print(ifile, '-dpng','-r300')
    crop_image(ifile)
    
    %%%%%%%%%%%
    figure(3)
    clf
    
    s = ["Bransfield" "Elephant" "Joinville" "West"];
    plot_standard_map(strata, 'centrePoint', [-58 -62], 'radius', 4, ...
        'strata', s, 'showStrataNames', false, ...
        'coastDetail', 'high')
    maxNASC = max(nasc.NASC);
    
    h = [];
    for i = 1:length(s)
        j = find(nasc.Stratum == s(i));
        h(i) = m_scatter(nasc.Longitude(j), nasc.Latitude(j), nasc.NASC(j)/maxNASC*200+1, 'filled', 'o');
    end
    
    m_grid('box', 'on')
    legend(h, s, 'Location', 'SouthEast')
    
    ifile = fullfile(resultsDir, 'NASC - AMLR.png');
    print(ifile, '-dpng','-r300')
    crop_image(ifile)
    
    %%%%%%%%%%%
    figure(4)
    clf
    
    s = ["ESS" "Sand" "SG" "SS" "AP" "SSI"];
    plot_standard_map(strata, 'centrePoint', [-45 -60], 'radius', 17.5, ...
        'strata', s, 'showStrataNames', false, ...
        'coastDetail', 'intermediate')
    maxNASC = max(nasc.NASC);
    
    h = [];
    for i = 1:length(s)
        j = find(nasc.Stratum == s(i));
        h(i) = m_scatter(nasc.Longitude(j), nasc.Latitude(j), nasc.NASC(j)/maxNASC*200+1, 'filled', 'o');
    end
    
    m_grid('box', 'on')
    legend(h, s, 'Location', 'SouthEast')
    
    ifile = fullfile(resultsDir, 'NASC - CCAMLR 2000.png');
    print(ifile, '-dpng','-r300')
    crop_image(ifile)
    
    %%%%%%%%%%%
    figure(5)
    clf
    
    s = ["SOI" "SOC" "SOF"];
    plot_standard_map(strata, 'centrePoint', [-45 -61], 'radius', 2.5, ...
        'strata', s, 'showStrataNames', false, ...
        'coastDetail', 'fine')
    maxNASC = max(nasc.NASC);
    
    h = [];
    for i = 1:length(s)
        j = find(nasc.Stratum == s(i));
        h(i) = m_scatter(nasc.Longitude(j), nasc.Latitude(j), nasc.NASC(j)/maxNASC*200+1, 'filled', 'o');
    end
    
    m_grid('box', 'on')
    legend(h, s, 'Location', 'SouthEast')
    
    ifile = fullfile(resultsDir, 'NASC - South Orkney.png');
    print(ifile, '-dpng','-r300')
    crop_image(ifile)
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Map of transect names. Bit complicated with lots of subplots
    figure(6)
    clf
    
    % Wide area strata
    subplot1(2,2, 'Gap', [0.05 0.05])
    subplot1(1)
    s = ["SS" "AP" "ESS"];
    plot_standard_map(strata, 'strata', s, 'showStrataNames', true, ...
        'coastDetail', 'low')
    plot_transect_names(nasc, s)
    m_grid('box', 'on')
   
    % South Shetland Islands
    %     subplot(2,2,2)
    %     s = ["SSI"];
    %     plot_standard_map(strata, 'strata', s, 'showStrataNames', true, ...
    %         'coastDetail', 'low', 'centrePoint', [-58.2 -62], ...
    %         'radius', 4)
    %     plot_transect_names(nasc, s)
    %     m_grid('box', 'on')

    % AMLR transects
     %    subplot(3,3,6)
         s = ["West" "Elephant" "Joinville" "Bransfield"];
         plot_standard_map(strata, 'strata', s, 'showStrataNames', true, ...
             'coastDetail', 'low', 'centrePoint', [-58.2 -62], ...
             'radius', 4)
         plot_transect_names(nasc, s)
         m_grid('box', 'on')
    
    % South Orkney #1
    
    subplot1(2)
    s = ["SOF"];
    plot_standard_map(strata, 'strata', s, 'showStrataNames', true, ...
        'coastDetail', 'high', 'centrePoint', [-45.5 -60.5], ...
        'radius', 2.5)
    plot_transect_names(nasc, s)
    m_grid('box', 'on')
    
    % South Orkney #2
    subplot1(3)
    s = ["SOC"];
    plot_standard_map(strata, 'strata', s, 'showStrataNames', true, ...
        'coastDetail', 'fine', 'centrePoint', [-46.5 -60.25], ...
        'radius', .8)
    plot_transect_names(nasc, s)
    m_grid('box', 'on')
    
    % South Georgia
    subplot1(4)
    s = ["SG"];
    plot_standard_map(strata, 'strata', s, 'showStrataNames', true, ...
        'coastDetail', 'high', 'centrePoint', [-36.3 -54], ...
        'radius', 2)
    plot_transect_names(nasc, s)
    m_grid('box', 'on')
    
    % South Sandwich Islands
%     subplot(3,3,7)
%     s = ["ESS" "Sand"];
%     plot_standard_map(strata, 'strata', s, 'showStrataNames', true, ...
%         'coastDetail', 'low', 'centrePoint', [-26 -57], ...
%         'radius', 5)
%     plot_transect_names(nasc, s)
%     m_grid('box', 'on')
    
    ifile = fullfile(resultsDir, 'Transects - as done and labelled.png');
    print(ifile, '-dpng','-r300')
    crop_image(ifile)
end

