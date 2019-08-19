%%%
% Apply the 3-requency dB-difference method to data from Kronprins Haakon
% and save in a form for the do_estimate_biomass.m script to use.

do_define_directories
dataDir = fullfile(baseDir, 'data', 'echo-integration', 'KPH-dBdiff');
dataDirTransects = fullfile(baseDir, 'data', 'echo-integration', 'KPH');

clear nasc
for i = ["038", "120", "200"]
    disp("Loading " + i + " kHz data")
    d = dir(fullfile(dataDir, 'ListUserFile05__F' + i + '*.txt'));
    data = read_lsss_output(5, fullfile(d(1).folder, d(1).name), true);
    % trim to 250 m max depth
    maxDepth = 250; % [m]
    numRows = maxDepth / data.pelagic_layer_thickness(1); % assume that thickness is the same for all cells
    data.pelagic_sa = data.pelagic_sa(1:numRows, :);
    
    nasc.("f" + i) = data.pelagic_sa';
    nasc.Latitude = data.lat';
    nasc.Longitude = data.lon';
    nasc.Ping_timestamp = datetime(data.timestamp', 'ConvertFrom', 'datenum');
end

% and filter the nasc based on transect on/off times and also apply
% transect/strata labels

transects = readtable(fullfile(dataDirTransects, 'transects.csv'));
% transects.csv files have the time with milliseconds and
% timezone and readtable returns the time as a character
% string, so deal with that (and truncates the milliseconds)
transects.Time = datetime(transects.Time, 'InputFormat', 'HH:mm:ss.SSSZ', 'TimeZone', 'UTC');
transects.Time.TimeZone = '';
transects.Time.Second = floor(transects.Time.Second);

% Combine transect date and time into a timestamp
transects.timestamp = transects.Date + timeofday(transects.Time);

% find the on/off times
tt = find(startsWith(transects.Event, "On"));

onTransect = [];
transectNames = [];
stratumNames = [];

% for each on/off period, work out the stratum and transect name and
% then accumulate this all into one structure
for i = 1:length(tt)
    onTime = transects.timestamp(tt(i));
    offTime = transects.timestamp(tt(i)+1);
    event = transects.Event(tt(i));
    
    % work out the stratum and transect
    t = split(event, ' ');
    [stratum, transect] = sortOutNaming(string(t{3}));
            
    % Find data within the current time period
    j = find(nasc.Ping_timestamp >= onTime & nasc.Ping_timestamp <= offTime);
    onTransect = [onTransect; j];
    
    % Make up transect and strata for each nasc value in the current
    % time period
    transectNames = [transectNames; repmat(transect, length(j), 1)];
    stratumNames = [stratumNames; repmat(stratum, length(j), 1)];
end

% and now sort out the nasc structure
nasc.Latitude = nasc.Latitude(onTransect);
nasc.Longitude = nasc.Longitude(onTransect);
nasc.Ping_timestamp = nasc.Ping_timestamp(onTransect);
nasc.Transect = transectNames;
nasc.Stratum = stratumNames;
nasc.Vessel = repmat("KPH", size(transectNames, 1), 1);
nasc.f038 = nasc.f038(onTransect,:);
nasc.f120 = nasc.f120(onTransect,:);
nasc.f200 = nasc.f200(onTransect,:);


% Calculate the dB difference windows. Results are almost identical to
% Table 3 in EMM-16/38. Small differences appear to be due to the rounding in
% that table (to 1 decimal place) and also that the TS values in the
% sdwba_ts.csv file (which come from Table 2 in the same document) are only
% to 2 decimal places. 
%ts = readtable(fullfile(resultsDir, 'sdwba_ts.csv'));
load(fullfile(resultsDir, 'ASAM_TS_krill_length_values_alt_fin'), 'T_TS', 'krill_ls')
ts.ts038 = T_TS(1,:)';
ts.ts120 = T_TS(2,:)';
ts.ts200 = T_TS(3,:)';
ts.length = krill_ls';

k = 1;
for len_min = (10:10:50)
    for len_max = len_min+10:10:60
        sizes = len_min:len_max;
        dBdiff1 = nan(size(sizes));
        dBdiff2 = nan(size(sizes));
        i = 1;
        for l = sizes
            min_i = find(ts.length == l);
            dBdiff1(i) = ts.ts120(min_i) - ts.ts038(min_i);
            dBdiff2(i) = ts.ts200(min_i) - ts.ts120(min_i);
            i = i + 1;
        end
        range1 = [min(dBdiff1) max(dBdiff1)];
        range2 = [min(dBdiff2) max(dBdiff2)];
        dbWindows.minLen(k) = sizes(1);
        dbWindows.maxLen(k) = sizes(end);
        dbWindows.range120_038(k,:) = range1;
        dbWindows.range200_120(k,:) = range2;
        k = k + 1;
    end
end

% Get lf's and work out what dB differences we are to use
load(fullfile(resultsDir, 'Trawls - data'), 'lf_raw', 'lf');
% Use the full set of trawls, not just those done by KPH. This makes things
% more comparable between analyses.

% Work out the dB difference windows for both the per-stratum and clusters length distributions.
for field = ["cluster", "strata"]
    for i = 1:length(lf.(field))
        toRemove = round(0.025 * lf.(field)(i).numLengths); % from each end
        trimmedLengths = sort(lf.(field)(i).lengths);
        trimmedLengths = trimmedLengths(toRemove+1:end-toRemove-1);
        [N, edges] = histcounts(trimmedLengths, 10:10:70);
        l = find(N > 0);
        lf.(field)(i).dBdiffLengthRange = [edges(l(1)) edges(l(end))];
        j = find(dbWindows.minLen == edges(l(1)) & ...
            dbWindows.maxLen == edges(l(end)));
        if length(j) ~= 1
            error('do_db_difference, when finding db Window')
        end
        lf.(field)(i).dbWindow120_038 = dbWindows.range120_038(j,:);
        lf.(field)(i).dbWindow200_120 = dbWindows.range200_120(j,:);
    end
end

% Filter the nasc using the dB differences
% do the per-stratum db differences
nasc.NASC_per_stratum = nan(size(nasc.f120));
%nasc.NASC_per_cluster = nan(size(nasc.f120));

% the dB differences
diff120_038 = 10*log10(nasc.f120)-10*log10(nasc.f038);
diff200_120 = 10*log10(nasc.f200)-10*log10(nasc.f120);

s = unique(nasc.Stratum);
for i = 1:length(s)
    % find the right dB window for this stratum
    s_i = find(strcmp({lf.strata.stratum}, s(i)));
    if ~isempty(s_i) % some strata are not in this dataset or aren't strata that we do biomass on 
        w1 = lf.strata(s_i).dbWindow120_038;
        w2 = lf.strata(s_i).dbWindow200_120;
    
        jj = nasc.Stratum == s(i) & ...
            diff120_038 >= w1(1) & diff120_038 <= w1(2) & ...
            diff200_120 >= w2(1) & diff200_120 <= w2(2);
        
        nasc.NASC_per_stratum(jj) = nasc.f120(jj);
    end
end

% and then collapse the depth axis and resample into 1 n.mi horizontal
% bins (as per CCAMLR procedure) - it currently in  m bins
nasc.NASC_per_stratum = nansum(nasc.NASC_per_stratum, 2);


% Format for the do_estimate_biomass script
nasc = rmfield(nasc, {'f038' 'f120' 'f200'});

nasc = struct2table(nasc);

% flag each nasc value to show if it was taken within (or not) civil
% daylight hours.
for i = 1:height(nasc)
    nasc.civilDaytime(i) = dayOrNight(nasc.Longitude(i), nasc.Latitude(i), nasc.Ping_timestamp(i));
end

% and save the combined dataset
save(fullfile(resultsDir, 'NASC - data - dB difference - KPH'), 'nasc')

% and do some comparison plots of NASC
% load(fullfile(resultsDir, 'NASC - data - dB difference - KPH'), 'nasc')
% nasc_diff = nasc;
% load(fullfile(resultsDir, 'NASC - data'), 'nasc')
% 
% s = unique(nasc_diff.Stratum+nasc_diff.Transect);

% figure(1)
% clf
% for i = 1:length(s)
%     i_diff = find(nasc_diff.Stratum + nasc_diff.Transect == s(i) & ...
%     nasc_diff.Vessel == "KPH");
%     i_swarm = find(nasc.Stratum + nasc.Transect == s(i) & nasc.Vessel == "KPH");
%     plot(nasc.Latitude(i_swarm), nasc.NASC(i_swarm), '.-')
%     hold on
%     plot(nasc_diff.Latitude(i_diff), nasc_diff.NASC_per_stratum(i_diff), '.-')
%     pause
%     clf
% end




