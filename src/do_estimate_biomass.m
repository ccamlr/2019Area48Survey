%%
% Estimate the krill biomass per strata and per survey

baseDir = 'I:\KRILL2019';
repoDir = '2019Area48SurveyRepo';
resultsDir = fullfile(baseDir, repoDir, 'results');

% Strata areas. 

% ** Should be calculated directly from the latest strata boundaries **

% CCAMLR 2000 strata areas from Bo workshop 2004, section 2.3
% AMLR strata areas from Figure 1 of EMM-11/26. 
% Units are km^2
s = {'AP', 473318; 'SS' 1109789; 'ESS', 321800; ...
    'SSI', 48654; 'SOI', 24409; 'SG', 25000; 'Sand', 62274; ...
    'Joinville', 18151; 'Elephant', 43865; 'West', 38524; 'South', 24479};
% and turn this into a more useful structure
clear strata_area
for i = 1:size(s,1)
    strata_area(i).name = s{i,1};
    strata_area(i).area = s{i,2};
end

% Length to weight relationship. No data currently available from the 2019
% survey, so use here the relationship for the 2000 survey (as per
% EMM-16/38, equation 7.

w = @(x) 2.236e-6 * (x*1e3).^3.314; % takes lengths in m and returns weight in grams.


% Krill TS values
load(fullfile(resultsDir, 'SDWBA-TS-2019'), 'krill_ts')
krill_len = [krill_ts.ts(:).ActualLength]; % [m]
krill_sigma = [krill_ts.ts(:).sigma_avg] *4*pi; % Why oh why are the CCAMLR equations not using backscattered sigma????

% NASC data.
load(fullfile(resultsDir, 'NASC - data'), 'nasc')

% Lengths per strata/cluster
load(fullfile(resultsDir, 'Trawls - data'), 'lf')

% Weighted conversion factor for each strata
for i = 1:length(lf.strata)
    
    histedges = lf.strata(i).histedges(1:end-1) * 1e-3; % [m]
    j = ismember(krill_len, histedges);
    sigma = krill_sigma(j); % [m^2]
    
    % Weighted conversion factor as per EMM-16/38, eqn 8. 
    C = sum(lf.strata(i).histcounts .* w(histedges)) ./ ...
         sum(lf.strata(i).histcounts .* sigma); % [g/m^2]
    lf.strata(i).C = C;
end

% Weighted conversion factor for each cluster
for i = 1:length(lf.cluster)
    
    histedges = lf.cluster(i).histedges(1:end-1) * 1e-3; % [m]
    j = ismember(krill_len, histedges);
    sigma = krill_sigma(j);
    
    % Weighted conversion factor as per EMM-16/38, eqn 8.
    C = sum(lf.cluster(i).histcounts .* w(histedges)) ./ ...
         sum(lf.cluster(i).histcounts .* sigma); % [g/m^2]
    lf.cluster(i).C = C;
end

% Associate conversion factor with each NASC value per strata
strata = {lf.strata.stratum};

% This is the slow way, but will catch nasc rows with unknown strata...
for i = 1:height(nasc)
    j = find(strcmp(nasc.Stratum(i), strata));
    if ~isempty(j)
        nasc.C(i) = lf.strata(j).C;
    else
        nasc.C(i) = NaN; % for nasc values in strata without a krill lf, or where a strata is not present.
    end
end

% Now calculate biomass density for each transect, stratum and finally
% survey. This follows the equations in EMM-16/38 section 9.

clear biomass
for k = 1:length(strata) % loop over all strata for which we have 
    biomass.strata(k).name = strata(k);
    transects = unique(nasc.Transect(nasc.Stratum == strata(k)));
    for j = 1:length(transects) % loop over all transects in each strata
        intervals = find(nasc.Transect == transects(j) & nasc.Stratum == strata(k));
        
        [rho_j, L_j] = calcDensity(nasc(intervals,:));
                
        biomass.strata(k).transect(j).name = transects(j);
        biomass.strata(k).transect(j).krillDensity = rho_j; % [g/m^2]
        biomass.strata(k).transect(j).length = L_j; % [nmi]
    end
    
    % Per stratum calculations
    
    % extract the per transect data from the biomass structure for
    % convenience.
    L_j = [biomass.strata(k).transect.length]; % [nmi]
    rho_j = [biomass.strata(k).transect.krillDensity]; % [g/m^2]
    N_k = length(transects); % [1]
    a = find(strcmp({strata_area.name}, strata(k)));
    A_k = strata_area(a).area; % [km^2]
    
    % Compute normalized transect weighting factors
    w_j = L_j / mean(L_j); % equation 11
        
    % Compute mean krill biomass density for entire survey area
    %rho_k = sum(w_j .* rho_j) / sum(w_j); % as in calcDensity script
    rho_k = sum(w_j .* rho_j) / N_k; % equation 14 [g/m^2]

    % Variance component
    VarComp_k = w_j.^2 .* (rho_j - rho_k).^2; % equation 13 

    % Compute variance of mean krill biomass density
    % Equation 15
    var_rho_k = N_k / (N_k-1) ...
        * sum(w_j.^2 .* (rho_k-rho_j).^2) / (sum(w_j)).^2 ...
        * sum(w_j.^2 .* (rho_k-rho_j).^2) / (N_k * (N_k-1));

    CV_k = 100 * sqrt(var_rho_k) / rho_k; % equation 16
    
    % and store the above values in a nicer form
    biomass.strata(k).area = A_k;
    biomass.strata(k).stratumLength = sum(L_j);
    biomass.strata(k).meanDensity = rho_k;
    biomass.strata(k).densityVariance = var_rho_k;
    biomass.strata(k).varianceComponent = VarComp_k;
    biomass.strata(k).CV = CV_k;
end
% and then calculate results for the entire survey
% this is a little tricky because there are multiple surveys (and
% overlapping strata), so we define which strata combine for each 'survey'.

biomass.survey(1).name = 'CCAMLR 2000';
biomass.survey(1).strata = ["AP", "SS", "SSI", "SOI", "Sand", "ESS"];

biomass.survey(2).name = 'AMLR';
biomass.survey(2).strata = ["Joinville", "Elephant", "West", "South"];

for s = 1:length(biomass.survey)
    % identify which strata are in this survey
    k = find(ismember([biomass.strata.name],biomass.survey(s).strata));
    
    % extract some per strata numbers into variables that match those used
    % in the equations
    A_k = [biomass.strata(k).area]; % [km^2]
    rho_k = [biomass.strata(k).meanDensity]; % [g/m^2]
    var_rho_k = [biomass.strata(k).densityVariance];
    VarComp_k = [biomass.strata(k).varianceComponent];
    
    % overall mean areal krill biomass density.
    rho = sum(A_k.*rho_k) / sum(A_k); % equation 17 [g/m^2]
    
    % overall survey variance
    var_rho = sum(A_k.^2 .* var_rho_k.^2) / sum(A_k.^2); % equation 18
    
    % overall CV of mean areal krill biomass density
    CV_rho = 100 * sqrt(var_rho) / rho; % equation 19
    
    % Biomass. 
    B_0 = sum(A_k .* rho_k); % equation 20 [tonnes]
    
    % Overall survey variance of B_0
    var_B_0 = sum(VarComp_k); % equation 21
    
    % Overall coefficient of variation of B_0
    CV_B_0 = 100 * sqrt(var_B_0) / B_0; % equation 22
    
    % and store the results in a nicer form
    biomass.survey(s).meanDensity = rho;
    biomass.survey(s).meanDensityVariance = var_rho;
    biomass.survey(s).meanDensityCV = CV_rho;

    biomass.survey(s).biomass = B_0;
    biomass.survey(s).variance = var_B_0;
    biomass.survey(s).CV = CV_B_0;
end


function [rho_j, L_j] = calcDensity(data)
   % calculate transect density.

   % Find extent of transect
   [~, i1] = max(data.Latitude);

   [~, i2] = min(data.Latitude);
   
   t = struct('lat1', data.Latitude(i1), 'lon1', data.Longitude(i1), ...
       'lat2', data.Latitude(i2), 'lon2', data.Longitude(i2));

   % Find the absolute difference between the start and end latitude of the transect
   delta_initial = abs(t.lat1 - t.lat2);
    
   % Find the distance in nautical miles between start and end points using
   % the Haversine formula
   [~, nmi] = haversine([t.lat1 t.lon1], [t.lat2 t.lon2]);
    
   % Compute expected change in latitude for each nautical mile
   e_d_lat = delta_initial/nmi ;
    
   % Total # of intervals
   numInts = height(data);
    
   W_I = nan(numInts, 1); % Interval lengths (n.mi.)
    
   W_I = ones(numInts, 1); %%%%%%%%%%%%%%% Do this properly!!!
    
   % If deviation from standard track line is < 10%, then simply set
   % weighting factor to 1
   W_I(W_I >= 0.9) = 1;
    
   % Compute transect length (n.mi.) by summing interval lengths
   L_j = sum(W_I);
    
   % Convert NASC to krill biomass density using weighting and conversion
   % factors. 
   % This is equation 9 in EMM-16/38 (but lacks the 1852^2 multplier).
   rho_j = sum(data.NASC / 1852^2.* data.C .* W_I) ./ L_j;
end


function [km, nmi, mi] = haversine(loc1, loc2)
    % HAVERSINE Compute distance between locations using Haversine formula
    % KM = HAVERSINE(LOC1, LOC2) returns the distance KM in km between
    % locations LOC1 and LOC2 using the Haversine formula. LOC1 and LOC2 are
    % latitude and longitude coordinates that can be expressed as decimal
    % degrees (where negative indicates West/South).
    %
    % [KM, NMI, MI] = HAVERSINE(LOC1, LOC2) returns the computed distance in
    % kilometers (KM), nautical miles (NMI), and miles (MI).
    %
    %
    % The first element indicates the latitude while the second is the
    % longitude.
    %
    % Notes
    % The Haversine formula is used to calculate the great-circle
    % distance between two points, which is the shortest distance over
    % the earth's surface.
    %
    % This program was created using equations found on the website
    % http://www.movable-type.co.uk/scripts/latlong.html
    % Created by Josiah Renfree
    % May 27, 2010
    %% Check user inputs
    % If two inputs are given, display error
    if ~isequal(nargin, 2)
        error('User must supply two location inputs')
        % If two inputs are given, handle data
    else
        locs = {loc1 loc2}; % Combine inputs to make checking easier
    end
    
    % Check that both cells are a 2-valued array
    if any(cellfun(@(x) ~isequal(length(x),2), locs))
        error('Incorrect number of input coordinates')
    end
    
    % Convert all decimal degrees to radians
    locs = cellfun(@(x) x .* pi./180, locs, 'UniformOutput', 0);
    
    % Begin calculation
    R = 6371; % Earth's radius in km
    delta_lat = locs{2}(1) - locs{1}(1); % difference in latitude
    delta_lon = locs{2}(2) - locs{1}(2); % difference in longitude
    a = sin(delta_lat/2)^2 + cos(locs{1}(1)) * cos(locs{2}(1)) * ...
        sin(delta_lon/2)^2;
    c = 2 * atan2(sqrt(a), sqrt(1-a));
    km = R * c; % distance in km
    % Convert result to nautical miles and miles
    nmi = km * 0.539956803; % nautical miles
    mi = km * 0.621371192; % miles
end