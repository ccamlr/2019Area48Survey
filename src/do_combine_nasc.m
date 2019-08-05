%%
% A script to combine the Echoview integrals from each vessel into one
% dataset, with each NASC value labelled with vessel and transect name.

baseDir = 'I:\KRILL2019';
repoDir = '2019Area48SurveyRepo';

dataDir = fullfile(baseDir, 'data', 'echo-integration');
resultsDir = fullfile(baseDir, repoDir, 'results');

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
            % transect numbers (see EMM-11/26)
            transects = {
                '20190205T102300', '20190205T153100', "Elephant", "02";
                '20190205T172300', '20190206T083050', "Elephant", "01";
                '20190206T083726', '20190206T164000', "Joinville", "01";
                '20190206T192400', '20190207T020500', "Joinville", "02";
                '20190207T052600', '20190207T101316', "Joinville", "03";
                '20190207T110955', '20190207T162900', "Bransfield", "01";
                '20190207T194500', '20190207T234100', "Bransfield", "02";
                '20190208T025800', '20190208T061900', "Bransfield", "03";
                '20190208T101200', '20190208T133300', "Bransfield", "04";
                '20190208T161900', '20190208T203000', "Bransfield", "05";
                '20190208T231100', '20190209T021800', "Bransfield", "06";
                '20190209T050200', '20190209T092301', "Bransfield", "07";
                '20190209T113440', '20190209T153450', "West", "07";
                '20190209T180130', '20190209T221645', "West", "06";
                '20190210T003000', '20190210T054150', "West", "05";
                '20190210T081200', '20190210T135850', "West", "04";
                '20190210T161030', '20190210T211230', "West", "03"};
                
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
            t = string(t{2});
            
            if startsWith(t, "AP")
                n.Stratum = repmat("AP", size(n.NASC, 1), 1); % Antarctic Pennisula
                n.Transect = repmat(extractAfter(t,2), size(n.NASC, 1), 1);
            elseif startsWith(t, "ELE") % AMLR stratum
                n.Stratum = repmat("Elephant", size(n.NASC, 1), 1);
                n.Transect = repmat(extractAfter(t,3), size(n.NASC, 1), 1);
            elseif startsWith(t, "WEST")  % AMLR stratum
                n.Stratum = repmat("West", size(n.NASC, 1), 1);
                n.Transect = repmat(extractAfter(t,4), size(n.NASC, 1), 1);
            elseif startsWith(t, "SSI")
                n.Stratum = repmat("SSI", size(n.NASC, 1), 1);
                n.Transect = repmat(extractAfter(t,3), size(n.NASC, 1), 1);
            elseif startsWith(t, "SSA") || startsWith(t, "SSB") || startsWith(t, "SSC")
                n.Stratum = repmat("ESS", size(n.NASC, 1), 1);
                n.Transect = repmat(extractAfter(t,2), size(n.NASC, 1), 1);
            elseif startsWith(t, "SS") % Scotia Sea. Must come after SSI
                n.Stratum = repmat("SS", size(n.NASC, 1), 1);
                n.Transect = repmat(extractAfter(t,2), size(n.NASC, 1), 1);
            elseif startsWith(t, "SG")
                n.Stratum = repmat("SG", size(n.NASC, 1), 1);
                n.Transect = repmat(extractAfter(t,2), size(n.NASC, 1), 1);
            elseif startsWith(t, "SOF") % South Orkneys Fixed
                n.Stratum = repmat("SOF", size(n.NASC, 1), 1);
                n.Transect = repmat(extractAfter(t,3), size(n.NASC, 1), 1);
            elseif startsWith(t, "SOC") % South Orkney Concentrated
                n.Stratum = repmat("SOC", size(n.NASC, 1), 1);
                n.Transect = repmat(extractAfter(t,3), size(n.NASC, 1), 1);
            elseif startsWith(t, "SO") % South Okrney Islands. Must come after SOF and SOC
                n.Stratum = repmat("SOI", size(n.NASC, 1), 1);
                n.Transect = repmat(extractAfter(t,2), size(n.NASC, 1), 1);
            elseif startsWith(t, "SA483")
                n.Stratum = repmat("SA483", size(n.NASC, 1), 1);
                n.Transect = repmat(extractAfter(t,6), size(n.NASC, 1), 1);
            elseif startsWith(t, "Sand")
                n.Stratum = repmat("Sand", size(n.NASC, 1), 1);
                n.Transect = repmat(extractAfter(t,4), size(n.NASC, 1), 1);
            else
                n.Stratum = repmat("Unknown", size(n.NASC, 1), 1);
                n.Transect = repmat("Unknown", size(n.NASC, 1), 1);
            end
            
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

% and save the combined dataset
save(fullfile(resultsDir, 'NASC - data'), 'nasc')


%%%%%%%%%%%%%%%%
% Maps

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

print(fullfile(resultsDir, 'NASC - by vessel'), '-dpng','-r300')

% Map coloured by stratum
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

print(fullfile(resultsDir, 'NASC - by stratum'), '-dpng','-r300')


% Map of transect data and their names
figure(3)
clf
plot_standard_map(strata)

transects = unique(nasc.Stratum + nasc.Transect);

for i = 1:length(transects)
    k = find(nasc.Stratum+nasc.Transect == transects(i));
    m_plot(nasc.Longitude(k), nasc.Latitude(k), 'k.');
    hold on
    % and a label at the northernmost part
    [~, n_i] = max(nasc.Latitude(k));
    m_text(nasc.Longitude(k(n_i)), nasc.Latitude(k(n_i)), transects(i))
end

m_grid('box', 'on')

print(fullfile(resultsDir, 'Transects - as done and labelled'), '-dpng','-r300')


