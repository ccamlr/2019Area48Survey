%%
%
% Calculate a survey-wide sound speed and absorption value using all of the
% CTD's taken during the 2019 survey

baseDir = 'I:\KRILL2019';
repoDir = '2019Area48SurveyRepo';
resultsDir = fullfile(baseDir, repoDir, 'results');
dataDir = fullfile(baseDir, 'data', 'ctd');

% polygons that define the survey strata
strata = jsondecode(fileread(fullfile(baseDir, repoDir, 'map_data', 'survey strata.geojson')));

% Load in the CTD data from the various ships
ctd_raw = load_ctd_data(dataDir);

% Tidy up the data somewhat
ctd = ctd_raw;
for i = 1:length(ctd)
    disp(['Processing vessel: ' ctd(i).vessel ', index ' num2str(i) ])
    
    % remove data shallower than some value
    minDepth = 0.5; % [m]
    jj = find(ctd(i).d > minDepth);
    ctd(i).s = ctd(i).s(jj);
    ctd(i).d = ctd(i).d(jj);
    ctd(i).t = ctd(i).t(jj);

    % remove downcasts. Except for DY and KPH data, which is already done
    if ~(strcmp(ctd(i).vessel, 'DY') || strcmp(ctd(i).vessel, 'KPH'))
        [~, j] = max(ctd(i).d);
        ctd(i).s = ctd(i).s(j:end);
        ctd(i).d = ctd(i).d(j:end);
        ctd(i).t = ctd(i).t(j:end);
    end
    
    % calculate sound speed
    ctd(i).c = sw_svel(ctd(i).s, ctd(i).t, ctd(i).d); % [m/s]
    % and absorptions
    ctd(i).alpha038 = sw_absorption(38, ctd(i).s, ctd(i).t, ctd(i).d, 'fandg', 8); % [dB/km]
    ctd(i).alpha120 = sw_absorption(120, ctd(i).s, ctd(i).t, ctd(i).d, 'fandg', 8); % [dB/km]
    ctd(i).alpha200 = sw_absorption(200, ctd(i).s, ctd(i).t, ctd(i).d, 'fandg', 8); % [dB/km]

    % Demer 2004 give how they calculated an average c and absorption
    % for the survey area. We follow that.
    
    % 1. Calculate values at 10 m intervals (take a mean for each 10 m bin)
    binEdges = 0:10:250; % [m]
    
    [N, edges, bin] = histcounts(ctd(i).d, binEdges);
    ii = bin ~= 0; % will ignore bins with no data in them cause accumarray doesn't like it.
    binCentre = 0.5*(edges(1:end-1)+edges(2:end));
    ctd(i).subsampled_depth = binCentre(1:max(bin(ii)));
    
    ctd(i).subsampled_c = accumarray(bin(ii), ctd(i).c(ii), [], @mean)';
    ctd(i).subsampled_alpha038 = accumarray(bin(ii), ctd(i).alpha038(ii), [], @mean)';
    ctd(i).subsampled_alpha120 = accumarray(bin(ii), ctd(i).alpha120(ii), [], @mean)';
    ctd(i).subsampled_alpha200 = accumarray(bin(ii), ctd(i).alpha200(ii), [], @mean)';
    ctd(i).subsampled_s = accumarray(bin(ii), ctd(i).s(ii), [], @mean)';
    ctd(i).subsampled_t = accumarray(bin(ii), ctd(i).t(ii), [], @mean)';
    
    if length(ctd(i).subsampled_depth) ~= length(ctd(i).subsampled_c)
        % means that some depth bins had no data in them (cause the CTD
        % sample rate was slow and the trawl changed depth quickly). This
        % produces zeros in the subsampled data, so replace those with
        % interpolated values.
        jj = find(ctd(i).subsampled_c == 0);
        ctd(i).subsampled_c(jj) = interp1(ctd(i).d, ctd(i).c, ctd(i).subsampled_depth(jj), 'linear');
        ctd(i).subsampled_alpha038(jj) = interp1(ctd(i).d, ctd(i).alpha038, ctd(i).subsampled_depth(jj), 'linear');
        ctd(i).subsampled_alpha120(jj) = interp1(ctd(i).d, ctd(i).alpha120, ctd(i).subsampled_depth(jj), 'linear');
        ctd(i).subsampled_alpha200(jj) = interp1(ctd(i).d, ctd(i).calpha200, ctd(i).subsampled_depth(jj), 'linear');
        ctd(i).subsampled_s(jj) = interp1(ctd(i).d, ctd(i).s, ctd(i).subsampled_depth(jj), 'linear');
        ctd(i).subsampled_t(jj) = interp1(ctd(i).d, ctd(i).t, ctd(i).subsampled_depth(jj), 'linear');
    end
    
    % and calculate a mean value for each CTD, using the 1/r^2 weighting.
    w = 1 ./ [ctd(i).subsampled_depth].^2;
    ctd(i).mean_c = sum([ctd(i).subsampled_c] .* w) / sum(w);
    ctd(i).mean_alpha038 = sum([ctd(i).subsampled_alpha038] .* w) / sum(w);
    ctd(i).mean_alpha120 = sum([ctd(i).subsampled_alpha120] .* w) / sum(w);
    ctd(i).mean_alpha200 = sum([ctd(i).subsampled_alpha200] .* w) / sum(w);
end

% Calculate the weighted means

% 2. Take a depth weighted-mean as per equation (1) in Demer 2004. This
% is over all CTD's in the survey

% do the average in different regions
w = 1 ./ [ctd.subsampled_depth].^2;
    
clear r
r.mean_s = sum([ctd.subsampled_s] .* w) / sum(w);
r.mean_t = sum([ctd.subsampled_t] .* w) / sum(w);

r.mean_c = sum([ctd.subsampled_c] .* w) / sum(w);
r.mean_alpha038 = sum([ctd.subsampled_alpha038] .* w) / sum(w);
r.mean_alpha120 = sum([ctd.subsampled_alpha120] .* w) / sum(w);
r.mean_alpha200 = sum([ctd.subsampled_alpha200] .* w) / sum(w);

results = struct2table(r);

disp('Mean values over entire dataset')
disp(['  Temperature = ' num2str(r.mean_t, '%.1f'), ' degC'])
disp(['  Salinity = ' num2str(r.mean_s, '%.1f'), ' PSU'])
disp(['  Sound speed = ' num2str(r.mean_c, '%.f'), ' m/s'])
disp(['  Absorption at  38 kHz = ' num2str(r.mean_alpha038, '%.1f'), ' dB/km'])
disp(['  Absorption at 120 kHz = ' num2str(r.mean_alpha120, '%.1f'), ' dB/km'])
disp(['  Absorption at 200 kHz = ' num2str(r.mean_alpha200, '%.1f'), ' dB/km'])

% Calculate an average sound speed and absorption over the survey area that
% is independent of the spatial density of the CTD stations. Do this by
% interpolating the depth weighted mean sound speed and alpha120 from each
% CTD onto a 'uniform' grid (is lat/lon so not quite uniform in space) over
% the strata and then taking the mean of those sound speeds and alphas.
%
% This is slightly different from Demer because it is taking the average of
% an average. 

% Grid that covers the survey area
[X,Y] = meshgrid(-75:1:-15, -70:1:-50);
% NaN out areas outside the survey area in the uniform grid
in = false(size(X));
for i = 1:length(strata.features)
    poly = squeeze(strata.features(i).geometry.coordinates);
    in = in | inpolygon(X, Y, poly(:,1), poly(:,2));
end
X(~in) = NaN;
Y(~in) = NaN;
% Interpolate sound speed and alpha onto the grid
Fc = scatteredInterpolant([ctd.lon]', [ctd.lat]', [ctd.mean_c]', 'linear', 'none');
Zc = Fc(X,Y);
Falpha038 = scatteredInterpolant([ctd.lon]', [ctd.lat]', [ctd.mean_alpha038]', 'linear', 'none');
Zalpha038 = Falpha038(X,Y);
Falpha120 = scatteredInterpolant([ctd.lon]', [ctd.lat]', [ctd.mean_alpha120]', 'linear', 'none');
Zalpha120 = Falpha120(X,Y);
Falpha200 = scatteredInterpolant([ctd.lon]', [ctd.lat]', [ctd.mean_alpha200]', 'linear', 'none');
Zalpha200 = Falpha200(X,Y);
% and take the average
uniform.c = nanmean(Zc(:));
uniform.alpha038 = nanmean(Zalpha038(:));
uniform.alpha120 = nanmean(Zalpha120(:));
uniform.alpha200 = nanmean(Zalpha200(:));

disp('Mean values over interpolated dataset')
disp(['  Sound speed = ' num2str(uniform.c, '%.f'), ' m/s'])
disp(['  Absorption at  38 kHz = ' num2str(uniform.alpha038, '%.1f'), ' dB/km'])
disp(['  Absorption at 120 kHz = ' num2str(uniform.alpha120, '%.1f'), ' dB/km'])
disp(['  Absorption at 200 kHz = ' num2str(uniform.alpha200, '%.1f'), ' dB/km'])

% 3. Take a harmonic mean weighted by the krill depth distribution.
% leave this one for later...
    
% A plot of CTD positions and sound speed and alpha contours. Uses the m_map toolbox with the gshhs
% coastline.
clf

for k = 1:2
    subplot(1,2,k)

    plot_standard_map(strata, 'showStrataNames', false)
    m_plot([ctd.lon],[ctd.lat], 'k.')
        
    % A contour plot of sound speed over the survey area
    if k == 1
        F = scatteredInterpolant([ctd.lon]', [ctd.lat]', [ctd.mean_c]', 'linear', 'none');
        contourText = 'Sound speed [m/s]';
    else % and absorption
        F = scatteredInterpolant([ctd.lon]', [ctd.lat]', [ctd.mean_alpha120]', 'linear', 'none');
        contourText = 'Absorption at 120 kHz [dB/km]';
    end
    
    [X,Y] = meshgrid(-75:1:-15, -70:1:-50);
    Z = F(X,Y);
    [cs, h] = m_contour(X,Y,Z, 'LineWidth', 1);
    clabel(cs, h, 'fontsize', 9);
    
    m_grid('box', 'on')
    %title(contourText)
end

print(fullfile(resultsDir, 'CTD - station map'), '-dpng','-r300')

% save the results
save(fullfile(resultsDir, 'CTD - data'), 'ctd', 'results', 'uniform')


function ctd = load_ctd_data(dataDir)
% Function that loads in the CTD data. Each ship has it's own pecularities,
% so there is code specific to each ship in this function.
%
% Returns a structure array that has fields 'lat', 'lon', 's' (salinity,
% PSU), 't' (temperature, degC), 'd' (depth, m), 'timestamp' (Matlab
% datenums), and 'vessel'.

    j = 1;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Load in DY data
    disp('Loading RSS Discovery data')
    d = dir(fullfile(dataDir, 'RSS', '*.mat'));

    for i = 1:length(d)
        c = load(fullfile(d(i).folder, d(i).name));
        
        % some data has NaN's in it (often for the first depth bin and also
        % for all depths beyond the max depth that the CTD went to). Remove them.
        ii = ~(isnan(c.salin1) | isnan(c.temp1) | isnan(c.press));

        ctd(j) = struct('lat', c.lat, 'lon', c.lon, 's', c.salin1(ii), ...
            't', c.temp1(ii), 'd', c.press(ii), 'timestamp', datenum(c.gtime), ...
            'vessel', 'DY');
        j = j + 1;
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Load in Cabo de Hornos data
    disp('Loading Cabo de Hornos data')
    d = dir(fullfile(dataDir, 'CDH', '*.cnv'));

    % read in the trawl positions and start times
    pos = readtable(fullfile(dataDir, 'CDH', 'CTD station positions.csv'));
    pos.datenum = datenum(pos.timestamp,'yyyymmddTHHMM');

    js = j-1;
    
    for i = 1:length(d)
        c = ctd_rd(fullfile(d(i).folder, d(i).name), 'LAB');

        % each ctd file has many casts in it, so split it.
        %
        % It looks like each time the ctd is out of the water, the depth is
        % mostly negative, so use that to split
        betweenCastsDepth = 0.2; % [m]
        
        dp = c.depSM;
        dp(dp <= betweenCastsDepth) = NaN;
        start_i = find((diff([0;~isnan(dp(:))]) == 1) == 1);
        stop_i  = find((diff([0;~isnan(dp(:))]) == -1) == 1);
        
        % and ignore datasets that are only a few measurements
        k = find((stop_i - start_i) > 10);
        start_i = start_i(k);
        stop_i = stop_i(k);
        
        % split out the casts and store
        for k = 1:length(start_i)
            a = start_i(k);
            b = stop_i(k);
            
            % to get the position of each CTD, we assume that the order of
            % CTD's are chronological, as are the trawl positions. The
            % timestamps in the CTD and trawl data don't match up and the
            % difference changes with time. There are the same number of
            % CTD's and trawls, so assume that they match up...
            
            ctd(j) = struct('lat', pos.latitude(j-js), 'lon', pos.longitude(j-js), 's', c.sal00(a:b), ...
                't', c.tv290C(a:b), 'd', c.depSM(a:b), ...
                'timestamp', datenum(c.gtime) + c.timeH(a)/24, ...
                'vessel', 'CDH');
            j = j + 1;            
        end
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Load in Kronprins Haakon data
    disp('Loading Kronprins Haakon data')
    d = dir(fullfile(dataDir, 'KPH', '*.cnv'));
    for i = 1:length(d)
        c = ctd_rd(fullfile(d(i).folder, d(i).name), 'NMEA');

        ctd(j) = struct('lat', c.latitude, 'lon', c.longitude, 's', c.sal00, ...
            't', c.t068C, 'd', c.depSM, 'timestamp', datenum(c.gtime), ...
            'vessel', 'KPH');
        j = j + 1;
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Load in More Sodruzhestva data
    disp('Loading More Sodruzhestva data')
    % station positions obtained from the trawl data (CTD was on the trawl)
    % station_id, latitude, longitude
    pos = [1735	-62.09266667	-64.92333333
           1825	-64	-65.25
           1826	-62.76866667	-66.11866667
           1827	-61.452	-66.96683333
           1828	-60.0205	-67.82133333
           1736	-60.6025	-65.93433333
           1633	-61.19316667	-62.071
           1634	-62.30116667	-61.21733333];
       
    d = dir(fullfile(dataDir, 'MS', '*.cnv'));
    for i = 1:length(d)
        c = ctd_rd(fullfile(d(i).folder, d(i).name), 'LAB');
        % The station number is at the start of the filename
        stationId = str2double((d(i).name(1:end-4)));
        ii = find(stationId == pos(:,1));
        
        ctd(j) = struct('lat', pos(ii, 2), 'lon', pos(ii, 3), 's', c.sal00, ...
            't', c.tv290C, 'd', c.depSM, 'timestamp', datenum(c.gtime), ...
            'vessel', 'MS');
        j = j + 1;
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Load in Fu Rong Hai data
    disp('Loading Fu Rong Hai data')
    % Generated by the Matlab Import tool
    clear opts
    opts = delimitedTextImportOptions("NumVariables", 8);
    opts.DataLines = [1, Inf];
    opts.Delimiter = " ";
    opts.VariableNames = ["salinity", "temperature", "pressure", "day", "month", "year", "time", "unknown"];
    opts.VariableTypes = ["double", "double", "double", "double", "categorical", "double", "datetime", "double"];
    opts = setvaropts(opts, 7, "InputFormat", "HH:mm:ss");
    opts = setvaropts(opts, 5, "EmptyFieldRule", "auto");
    opts.ExtraColumnsRule = "ignore";
    opts.EmptyLineRule = "read";
    opts.ConsecutiveDelimitersRule = "join";
    opts.LeadingDelimitersRule = "ignore";

    % station (CTD/trawl) positions extracted from cruise log
    pos = readtable(fullfile(dataDir, 'FRH', 'CTD station positions.csv'));
    
    d = dir(fullfile(dataDir, 'FRH', '*.asc'));
    for i = 1:length(d)
        c = readtable(fullfile(d(i).folder, d(i).name), opts);
        
        % this vessel only provided data in Feb, so hard-code that...
        % and we only want one time per cast, so choose the time of the first measurement
        tt = datenum(c.year(1), 2, c.day(1), c.time.Hour(1), c.time.Minute(1), c.time.Second(1));
        
        % station name at the end of the filename
        station_id = d(i).name(end-9:end-4);
        ii = find(strcmp(station_id, pos.StationID) == 1);
        if length(ii) ~= 1
            error(['No station position found for CTD file ' fullfile(d(i).folder, d(i).name) ' - skipped'])
        else
            ctd(j) = struct('lat', pos.lat(ii), 'lon', pos.lon(ii), 's', c.salinity, ...
                't', c.temperature, 'd', c.pressure, 'timestamp', tt, ...
                'vessel', 'FRH');
            j = j + 1;
        end
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Load in Kwang Ja Ho data
    disp('Loading Kwang Ja Ho data')
    d = dir(fullfile(dataDir, 'KJH', '*.cnv'));

    % station (CTD/trawl) positions extracted from cruise log
    pos = readtable(fullfile(dataDir, 'KJH', 'CTD station positions.csv'));

    for i = 1:length(d)
        if strcmp(d(i).name(1:4), 'TNH1') || strcmp(d(i).name(1:4), 'CAL_') || strcmp(d(i).name(1:5), '05-07')
            % ignore these ones
        else
            c = ctd_rd(fullfile(d(i).folder, d(i).name), 'LAB');
            % station name is at the start of the filename
            station_id = d(i).name(1:5);
            ii = find(strcmp(station_id, pos.station) == 1);
            if length(ii) ~= 1
                warning(['No station position found for CTD file ' fullfile(d(i).folder, d(i).name) ' - skipped'])
                if strcmp(station_id(1:5), '05-02')
                    disp('  ->using position from cruise plan')
                    lon = -55.75;
                    lat = -60.25;
                elseif strcmp(station_id(1:5), '11-07')
                    disp('  ->using position from cruise plan')
                    lon = -58.5;
                    lat = -61.5;
                elseif strcmp(station_id(1:5), '12-08')
                    disp('  ->using position from cruise plan')
                    lon = -59.0;
                    lat = -61.75;
                else
                    
                end
            else
                lat = pos.lat(ii);
                lon = pos.lon(ii);
            end
            ctd(j) = struct('lat', lat, 'lon', lon, 's', c.sal00, ...
                't', c.tv290C, 'd', c.depSM, 'timestamp', datenum(c.gtime), ...
                'vessel', 'KJH');
            j = j + 1;
        end
    end
end

