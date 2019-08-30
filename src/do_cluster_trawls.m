%%
% Load in krill lengths from trawls and calculate length lf in various
% ways (by vessel, by strata, by cluster).

do_define_directories
dataDir = fullfile(baseDir, 'data', 'catch');
catchDir = fullfile(baseDir, repoDir, 'catch_data');

clear lf

% Load in the krill length data from the various ships
lf_raw = load_lf_data(dataDir);

% L50's from Bjørn Krafft, based on trawl type for each vessel
L50 = struct('vessel', ["KPH" "CDH" "FRH" "MS" "KJH" "RRS"], ...
    'L50', [15.01 15.01 31.92 41.55 25.72 12.96]);

% Some lf about the trawls per vessel
vessels = unique({lf_raw.vessel});
for i = 1:length(vessels)
    j = find(strcmp({lf_raw.vessel}, vessels{i}));
    ll = cat(1,lf_raw(j).lengths);
    
    k = find(L50.vessel == vessels(i));
    
    numTidStn = sum(strcmp({lf_raw(j).type}, 'TID'));
    
    lf.vessel(i) = struct('vessel', vessels{i}, ...
        'mean', mean(ll), ...
        'std', std(ll), 'min', min(ll), ...
        'max', max(ll), 'numLengths', length(ll), ...
        'numStations', length(j), ...
        'numTidStations', numTidStn, ...
        'allLengths', ll, 'L50', L50.L50(k));
end

% Some lf per station
lf.station.mean = arrayfun(@(x) mean(x.lengths), lf_raw);
lf.station.mean_text = cellstr(num2str(lf.station.mean', '%.f'))';
lf.station.std = arrayfun(@(x) std(x.lengths), lf_raw);
lf.station.num = arrayfun(@(x) length(x.lengths), lf_raw);

% Store the min and max lengths for each station.
for i = 1:length(lf.station.num)
    if lf.station.num(i) > 0
        lf.station.min(i) = min(lf_raw(i).lengths);
        lf.station.max(i) = max(lf_raw(i).lengths);
    else
        lf.station.min(i) = NaN;
        lf.station.max(i) = NaN;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Clustering, as per the EMM report (and hence similar to the 2000 CCAMLR method)
% Uses the Matlab statistics toolbox functions.

% only do this on trawls with more than minNumKrill krill
minNumKrill = 20;
stationsToUse = find(lf.station.num >= minNumKrill);
stationsToNotUse = find(lf.station.num < minNumKrill);

% update the vessel_lf structure with some more stats
for i = 1:length(lf.vessel)
    % stations for each vessel with more than 20 krill in them
    stns20 = lf.station.num >= minNumKrill & strcmp({lf_raw.vessel}, lf.vessel(i).vessel);
    lf.vessel(i).numClusteredStations = sum(stns20);
    % lf from these trawls
    lf.vessel(i).clusteredLengths = cat(1,lf_raw(stns20).lengths);
end

% combine all lf data into a length-binned single matrix.
min_l = min(lf.station.min(stationsToUse)); % min krill length
max_l = max(lf.station.max(stationsToUse)); % max krill length
lengths = min_l:1:max_l; % bins

X = zeros(length(lf_raw(stationsToUse)), length(lengths));
for i = 1:length(lf_raw(stationsToUse))
    bins = histcounts(lf_raw(stationsToUse(i)).lengths, [lengths lengths(end)+1]);
    if sum(bins) > 0
        X(i,:) = bins;
    else
        X(i,:) = nan(size(bins));
    end
end

% Make up a table that suits Tor Knutsen's R code for clustering. This iwll
% be written to a file later on.
tor = array2table(X');
% and sort out the row and column names
for i = 1:length(stationsToUse)
    sName = matlab.lang.makeValidName(num2str(lf_raw(stationsToUse(i)).station));
    tor.Properties.VariableNames{i} = [lf_raw(stationsToUse(i)).vessel sName];
end
for i = 1:length(lengths)
    tor.Properties.RowNames{i} = ['S' num2str(lengths(i))];
end

% Do the hierarchical clustering.

% log2 transform (log2(a+1) a = number of individuals in a length class,
% center on a zero mean, and standarise to unit variance
Xd = log2(X+1);
m = mean(Xd, 2);
s = std(Xd, 0, 2);
Xd = (Xd - repmat(m, 1, size(Xd, 2))) ./ repmat(s, 1, size(Xd, 2));

D = pdist(Xd, 'euclidean');

% Try several linkage method and pick the best
link_methods = {'average' 'centroid' 'complete' 'median' 'single' 'ward' 'weighted'};
aggloCoeff = zeros(size(link_methods));
for i = 1:length(link_methods)
    Z = linkage(D, link_methods{i});
    % how good was the clustering?
    aggloCoeff(i) = cophenet(Z, D);
    disp([link_methods{i} ' ' num2str(aggloCoeff(i))])
end
% Use the best one...
[~, i] = max(aggloCoeff);
Z = linkage(D, link_methods{i});
disp(['Best linkage method is: ' link_methods{i} ', coeff = ' num2str(aggloCoeff(i))])

clusters = cluster(Z, 'MaxClust',3);
% This does it in one go, but we don't get the data necessary to draw a
% dendrogram.
%T = clusterdata(X, 'Distance', 'euclidean', 'Linkage', 'ward', 'MaxClust', 3);


% lf per cluster
for i = 1:length(unique(clusters))
    j = clusters == i;
    ll = cat(1,lf_raw(stationsToUse(j)).lengths);
    N = histcounts(ll, lengths);
    lf.cluster(i) = struct('cluster', ['Cluster ' num2str(i)], ...
        'mean', mean(ll), 'std', std(ll), 'min', min(ll), ...
        'max', max(ll),  'numLengths', length(ll), ...
        'numStations', sum(j), 'stationIndex', stationsToUse(j), ...
        'histcounts', N, 'histedges', lengths, ...
        'lengths', ll);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% And now the lf per stratum

% Load survey strata
strata = jsondecode(fileread(fullfile(baseDir, repoDir, 'map_data', 'survey strata.geojson')));

% and calculate some lf
for i = 1:length(strata.features)
    poly = squeeze(strata.features(i).geometry.coordinates);
    in = inpolygon([lf_raw.lon], [lf_raw.lat], poly(:,1), poly(:,2));
    ll = cat(1, lf_raw(in).lengths);
    N = histcounts(ll, lengths);
    lf.strata(i) = struct('stratum', strata.features(i).properties.stratum, ...
        'mean', mean(ll), ...
        'std', std(ll), 'min', min(ll), ...
        'max', max(ll),  'numLengths', length(ll), ...
        'numStations', sum(in), ...
        'histcounts', N, 'histedges', lengths, ...
        'lengths', ll);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do the lf per the 3 main CCAMLR 2000 strata, as discussed in SG-ASAM-2019
% Load survey strata

% This was done by Keith Reid and is loaded here...
lff = readtable(fullfile(catchDir, 'jason_indeed.csv'));
lff = removevars(lff, {'Var1'});
Klen = [lff.Klen; lff.Klen(end)+1]; % HACK. Define the upper limit of the last bin as 1mm more

for i = 1:length(strata.features)
    st_name = strata.features(i).properties.stratum;
    if ismember(st_name, {'AP' 'SSI' 'West' 'Elephant' 'Bransfield' 'Joinville'}) % Using AP
        disp(['Stratum ' st_name ' given lf from AP'])
        lf.strata(i).ASAM2019_normalised_lf_prop = lff.AP';
        lf.strata(i).ASAM2019_normalised_lf_len = Klen';
    elseif ismember(st_name, {'SS' 'SG' 'SOI' 'SOC' 'SOF' 'WCB'}) % Using SS
        disp(['Stratum ' st_name ' given lf from SS'])
        lf.strata(i).ASAM2019_normalised_lf_prop = lff.SS';
        lf.strata(i).ASAM2019_normalised_lf_len = Klen';
    elseif ismember(st_name, {'ESS' 'Sand'}) % Using ESS
        disp(['Stratum ' st_name ' given lf from ESS'])
        lf.strata(i).ASAM2019_normalised_lf_prop = lff.ESS';
        lf.strata(i).ASAM2019_normalised_lf_len = Klen';
    else
        error('Unknown strata!')
    end
end

% His process is repeated here cause I want this processing to be in
% here...

% and we do a comparison to check.



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Output the results in various ways:

% Summary of trawls per vessel
disp('Trawls by vessel:')
disp(struct2table(lf.vessel))
disp(['Total number of trawls was ' num2str(sum([lf.vessel.numStations]))])
disp(['Total number of pre-determined trawls with >20 krill was ' num2str(sum([lf.vessel.numClusteredStations]))])
disp(['Total number of TID trawls was ' num2str(sum([lf.vessel.numTidStations]))])

% Summary of trawls per strata
disp('Trawls by strata:')
disp(struct2table(lf.strata))
% Summary of trawls per cluster
disp('Trawls by cluster:')
disp(struct2table(lf.cluster))

% Print out the overall mean length and std
ll = [];
for i = 1:length(lf_raw)
    ll = [ll; lf_raw(i).lengths];
end
disp(['Overall krill mean length: ' num2str(mean(ll)) ' mm'])
disp(['Overall krill std length: ' num2str(std(ll)) ' mm'])

% save the results from the processing
save(fullfile(resultsDir, 'Trawls - data'), 'lf_raw', 'lf', 'aggloCoeff');

% save the results in a form that Tor Knutsen can easily use for his
% clustering R code
writetable(tor, fullfile(resultsDir, 'Trawls - for Tor.csv'), ...
    'FileType', 'text', 'WriteRowNames', true);

% and save in a format for Keith Reid
lf_keith = lf_raw;

for i = 1:length(lf_raw)
    if isnan(lf_raw(i).timestamp)
        lf_keith(i).timestamp = 0;
    end
end
% then write out a flattened file

fid = fopen(fullfile(resultsDir, 'Trawls - for Keith.csv'), 'w');
fprintf(fid, 'vessel,station code,year,month,day,start time,lat,lon,length,total\n');
for i = 1:length(lf_keith)
    header = sprintf('%s,%s,%s,%f,%f', lf_keith(i).vessel, num2str(lf_keith(i).station), ...
        datestr(lf_keith(i).timestamp, 'yyyy,mm,dd,HHMM'), lf_keith(i).lat, lf_keith(i).lon);
    N = histcounts(lf_keith(i).lengths, lengths);
    for j = 1:length(N)
        fprintf(fid, '%s,%d,%d\n', header, lengths(j), N(j));
    end
end
fclose(fid);

% and save the lengths per stratum, for checking with Martin

fid = fopen(fullfile(resultsDir, 'Trawls - per strata.csv'), 'w');
fprintf(fid, 'stratum,lengths\n');
for i = 1:length(lf.strata)
    fprintf(fid, '%s', lf.strata(i).stratum);
    fprintf(fid, ',%d', lf.strata(i).lengths);
    fprintf(fid, '\n');
end

%%%%%%%%%%%%%%%%%%%%%%
% Do some plots
%%%%%%%%%%%%%%%%%%%%%%
figure(1) % Cluster lfs
clf
for i = 1:length(unique(clusters))
    subplot(length(unique(clusters)), 1, i)
    j = clusters == i;
    histogram(cat(1,lf_raw(stationsToUse(j)).lengths), lengths)
    textLoc(['Cluster ' num2str(i)], 'NorthWest');
end
xlabel('Length (mm)')

ifile = fullfile(resultsDir, 'Trawls - cluster lf.png');
print(ifile, '-dpng','-r300')
crop_image(ifile)

%%%%%%%%%%%%%%%%%%%%%%%
figure(2) % Dendrogram of clusters
clf
dendrogram(Z,size(Z,1))
set(gca,'XTickLabel', {}, 'YTickLabel', {})

ifile = fullfile(resultsDir, 'Trawls - dendrogram.png');
print(ifile, '-r300', '-dpng')
crop_image(ifile)

%%%%%%%%%%%%%%%%%%%%%%
figure(3) % Map of stations, clusters, and strata
clf
plot_standard_map(strata, 'showStrataNames', false)

% plot the station positions
clear h labels

% Plot the stations that weren't used
for i = 1:length(stationsToNotUse)
    h(1) = m_scatter([lf_raw(stationsToNotUse(i)).lon], ...
        [lf_raw(stationsToNotUse(i)).lat], 10, 'k');
end
labels{1} = 'Not used';
    
% Plot the stations that were used in the clustering
for i = 1:length(unique(clusters))
    j = clusters == i;
    ii = stationsToUse(j);
    h(i+1) = m_scatter([lf_raw(ii).lon], [lf_raw(ii).lat], 15, 'filled');
    labels{i+1} = ['C' num2str(i)];
end

m_grid('box', 'on')
h = legend(h, labels, 'Location', 'southeast');
title(h, 'Clusters')

ifile = fullfile(resultsDir, 'Trawls - map.png');
print(ifile, '-dpng','-r300')
crop_image(ifile)

%%%%%%%%%%%%%%%%
figure(4) % All the lf's
clf
for i = 1:length(lf.strata)
    subplot(5,4,i)
    histogram('BinEdges', lf.strata(i).histedges, ...
        'BinCounts', lf.strata(i).histcounts, 'EdgeColor', 'none', ...
        'FaceColor', 'k')
    textLoc(lf.strata(i).stratum, 'NorthWest');
    if i >= 14
        xlabel('Length (mm)')
    end
end
for i = 1:length(lf.cluster)
    j = length(lf.strata) + i;
    subplot(5,4,j)
    histogram('BinEdges', lf.cluster(i).histedges, ...
        'BinCounts', lf.cluster(i).histcounts, 'EdgeColor', 'none', ...
        'FaceColor', 'k')
    textLoc(lf.cluster(i).cluster, 'NorthWest');
    if j >= 14
        xlabel('Length (mm)')
    end
end

ifile = fullfile(resultsDir, 'Trawls - lf per stratum.png');
print(ifile, '-dpng', '-r300')
crop_image(ifile)

%%%%%%%%%%%%%
figure(5) % lf's per vessel
clf
bins = 10:67;
for i = 1:length(lf.vessel)
    subplot(3,2,i)

    histogram(lf.vessel(i).clusteredLengths, 'BinEdges', bins, ...
        'EdgeColor', 'none', ...
        'FaceColor', 'k')
    hold on
    plot([lf.vessel(i).L50 lf.vessel(i).L50], ...
        get(gca,'YLim'), 'k', 'LineWidth', 2)
    textLoc(lf.vessel(i).vessel, 'NorthEast');
    %textLoc(['N=' num2str(lf.vessel(i).numStations)], 'NorthEast');
    if i >= 5
        xlabel('Length (mm)')
    end
end

ifile = fullfile(resultsDir, 'Trawls - lf per vessel.png');
print(ifile, '-dpng', '-r300')
crop_image(ifile)


figure(6) % The lf's used for each strata, as agreed upon in the 2019 ASAM meeting
clf
for i = 1:length(lf.strata)
    subplot(4,4,i)
    histogram('BinEdges', lf.strata(i).ASAM2019_normalised_lf_len, ...
        'BinCounts', lf.strata(i).ASAM2019_normalised_lf_prop, 'EdgeColor', 'none', ...
        'FaceColor', 'k')
    textLoc(lf.strata(i).stratum, 'NorthWest');
    if i >= 10
        xlabel('Length (mm)')
    end
end

ifile = fullfile(resultsDir, 'Trawls - lf per stratum ASAM2019.png');
print(ifile, '-dpng', '-r300')
crop_image(ifile)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function lf = load_lf_data(dataDir)
    % Am only interested in the length-frequencies, so will reduce the
    % datasets to this if not already...

    lf = struct('vessel', [], 'station', [], 'lengths', [], ...
        'lat', [], 'lon', [], 'timestamp', [], 'type', []);
    
    % RRS Discovery
    disp('Loading RRS Discovery data')
    c = readtable(fullfile(dataDir, 'RRS.xlsx')); 
    s = readtable(fullfile(dataDir, 'RRS_station_info.csv'));
    
    stations = unique(s.EventNo);
    k = length(lf);
    for i = 1:length(stations)
        j = c.Event_number == stations(i);
        
        kk = find(s.EventNo == stations(i) & strcmpi(s.Action, 'Net 1 opened'));
        if length(kk) ~= 1
            warning('No single station info found')
        end
        
        if strcmp(s.Comment(kk), 'Stratified')
            stnType = 'STN';
        else
            stnType = 'TID';
        end
        
        lf(k) = struct('vessel', 'RRS', 'station', stations(i), ...
            'lengths', c.Length(j), ...
            'lat', s.Latitude(kk), 'lon', s.Longitude(kk), ...
            'timestamp', datenum(s.timestamp(kk)), ...
            'type', stnType);
        k = k + 1;
    end
    
    % More Sodruzhestva
    disp('Loading More Sodruzhestva data')
    c = readtable(fullfile(dataDir, 'MS.xlsx'), 'Sheet', 'look here gavin');
    s = readtable(fullfile(dataDir, 'MS.xlsx'), 'Sheet', 'station position');
    s.lat = (fix(s.lat/1e4) + (s.lat/1e4 - fix(s.lat/1e4))*100/60);
    s.lon = (fix(s.lon/1e4) + (s.lon/1e4 - fix(s.lon/1e4))*100/60);
    s.timestamp = datenum([num2str(s.date) num2str(s.time_UTC, '%04d00')], 'yyyymmddHHMMSS');
    
    k = length(lf) + 1;
    for st = 1:length(s.st_no)
        if s.st_no(st) == 1634 || s.st_no(st) == 1828
            lengths = c.(['st' num2str(s.st_no(st))]);
            lengths = lengths(~isnan(lengths));
        else
            lengths = [];
        end
    
        lf(k) = struct('vessel', 'MS', 'station', s.st_no(st), ...
            'lengths', lengths, ...
            'lat', s.lat(st), 'lon', s.lon(st), ...
            'timestamp', s.timestamp(st), ...
             'type', 'STN');
        k = k + 1;        
    end

    % Fu Rong Hai
    disp('Loading Fu Rong Hai data')
    c = readtable(fullfile(dataDir, 'FRH.xlsx')); % catch info
    s = readtable(fullfile(dataDir,'FRH_station_info.csv')); % station info
    
    stations = unique(c.Station);
    k = length(lf) + 1;
    for i = 1:length(stations) % loop over stations in catch info
        j = find(strcmp(c.Station, stations(i))); % all lines in catch info for current station
        
        kk = find(strcmp(c.Station(j(1)), s.station)); % station info for the current station
        if length(kk) ~= 1
           warning('No single station info found')
        end
        
        lf(k) = struct('vessel', 'FRH', 'station', c.SYNOPTIC_ST(j(1)), ...
            'lengths', c.Length(j), ...
            'lat', s.latitude(kk), 'lon', s.longitude(kk), ...
            'timestamp', datenum(s.timestamp(kk)), ...
             'type', 'STN');
        k = k + 1;
    end
    
    % Kronprins Haakon
    % Note. This dataset, obtained from the at-sea plankton database
    % differs a little from the set that Bjørn prepared for Gavin.
    disp('Loading Kronprins Haakon data') 
    c = readtable(fullfile(dataDir, 'KPH_krill_catch_extract.csv'));
    
    % fix some column data types
    c.lon = str2double(c.lon);
    c.lat = str2double(c.lat);
    c.len = str2double(c.len);
    
    stations = unique(c.serienr);
    k = length(lf) + 1;
    for i = 1:length(stations)
        j = find(c.serienr == stations(i));
        
        % pick out hour and minute of the time
        h = fix(c.starttid(j(1))/100);
        m = c.starttid(j(1)) - h*100;
                
        lengths = repelem(c.len(j), c.antall(j));
        lengths = reshape(lengths, [], 1); % force to be a column vector

        stnType = 'STN';
        if ismember(stations(i), [4002])
            stnType = 'TID';
        end
                
        if ~isempty(lengths)
            lf(k) = struct('vessel', 'KPH', 'station', stations(i), ...
                'lengths', lengths, ...
                'lat', c.lat(j(1)), 'lon', c.lon(j(1)), ...
                'timestamp', datenum(c.aar(j(1)), c.mnd(j(1)), c.dag(j(1)), ...
                h, m, 0),  'type', stnType);
            k = k + 1;
        end
    end
    
    % Cabo de Hornos
    % Note. This dataset, obtained from the at-sea plankton database,
    % differs a lot from the set that Bjørn prepared for Gavin.
    disp('Loading Cabo de Hornos data')
    c = readtable(fullfile(dataDir, 'CDH_krill_catch_extract.csv'));
    
    % fix some column data types
    c.lon = str2double(c.lon);
    c.lat = str2double(c.lat);
    c.len = str2double(c.len);
    
    stations = unique(c.serienr);
    k = length(lf) + 1;
    for i = 1:length(stations)
        j = find(c.serienr == stations(i));
        
        % pick out hour and minute of the time
        h = fix(c.starttid(j(1))/100);
        m = c.starttid(j(1)) - h*100;
        
        stnType = 'STN';
        if ismember(stations(i), [4266 4267 4273 4308 4311 4317])
            stnType = 'TID';
        end
        
        lengths = repelem(c.len(j), c.antall(j));
        
        if ~isempty(lengths)
            lf(k) = struct('vessel', 'CDH', 'station', stations(i), ...
                'lengths', lengths, ...
                'lat', c.lat(j(1)), 'lon', c.lon(j(1)), ...
                'timestamp', datenum(c.aar(j(1)), c.mnd(j(1)), c.dag(j(1)), ...
                h, m, 0),  'type', stnType);
            k = k + 1;
        end
    end
    
    % Kwang Ja Ho
    disp('Loading Kwang Ja Ho data')
    
    c = readtable(fullfile(dataDir, 'KJH.xlsx'), 'Sheet', 'Krill length frequency data');
    st = readtable(fullfile(dataDir, 'KJH.xlsx'), 'Sheet', 'Catch');
    s = readtable(fullfile(dataDir, 'KJH_station_info.csv'));
    
    % remove the entries that are thought to be ice krill, not Antarctic
    % krill
    i = ~strcmp(c.comments, 'ice');
    c = c(i,:);
    
    % correct some inconsistencies in the station names so they match those
    % in the s table.
    i = strcmp(c.station, '05.05-08(158)');
    c.station(i) = {'05.5-08(158)'};
    % and trim off the station number in the c table, leaving just the AMLR station name
    c.station = cellfun(@(S) S(1:end-5), c.station, 'UniformOutput', false);
    % trim off trailing spaces...
    c.station = strtrim(c.station);
    % the hyphen is not a hyphen in the data from the
    % spreadsheet (probably due to the Korean coding of the spreadsheet
    % file), so fix that.
    c.station = regexprep(c.station, '\xAD', '-');
    st.Station = regexprep(st.Station, '\xAD', '-');
    
    % make a datetime from the separate fields in the st table
    st.startTimestamp = datetime(num2str(st.Date), 'InputFormat', 'yyyyMMdd') + days(st.StartTime);
    st.stopTimestamp = datetime(num2str(st.Date), 'InputFormat', 'yyyyMMdd') + days(st.EndTime);
    
    stations = unique(c.station);
    k = length(lf) + 1;
    for i = 1:length(stations)
        j = strcmpi(c.station, stations(i));
        
        % Find the matching station in the s table (to get lat/lon)
        kk = find(strcmp(s.station, stations(i)));
        if length(kk) ~= 1
            warning(['No info found for station ' stations{i}])
        end
        
        % Find the matching station in the st table (to get start/stop
        % times)
        jj = strcmpi(st.Station, stations(i));

        lengths = c.length_mm_(j);
        
        if ~isempty(lengths)
            lf(k) = struct('vessel', 'KJH', 'station', stations(i), ...
                'lengths', lengths, ...
                'lat', s.lat(kk), 'lon', s.lon(kk), ...
                'timestamp', datenum(st.startTimestamp(jj)), ...
                'type', s.type(kk));
            k = k + 1;
        end
    end
    
end
