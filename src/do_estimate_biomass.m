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
krill_sigma = [krill_ts.ts(:).sigma_avg]*4*pi; % [m^2] Why oh why are the CCAMLR equations not using backscattered sigma????

% NASC data.
load(fullfile(resultsDir, 'NASC - data'), 'nasc')

% Lengths per strata/cluster
load(fullfile(resultsDir, 'Trawls - data'), 'lf')

% Weighted conversion factor for each strata and cluster
for t = ["strata" "cluster"]
    for i = 1:length(lf.(t))
        histedges = lf.(t)(i).histedges(1:end-1) * 1e-3; % [m]
        j = ismember(krill_len, histedges);
        sigma = krill_sigma(j); % [m^2]
        
        % Weighted conversion factor as per EMM-16/38, eqn 8.
        C = sum(lf.(t)(i).histcounts .* w(histedges)) ./ ...
            sum(lf.(t)(i).histcounts .* sigma); % [g/m^2]
        lf.(t)(i).C = C;
    end
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
% survey. This follows the equations in EMM-16/38 section 9, with some
% corrections.

clear biomass
for k = 1:length(strata) % loop over all strata 
    biomass.strata(k).name = strata(k);
    transects = unique(nasc.Transect(nasc.Stratum == strata(k)));
    for j = 1:length(transects) % loop over all transects in each strata
        % Find all NASC values for the current transect
        intervals = find(nasc.Transect == transects(j) & nasc.Stratum == strata(k));
        % Calculate the mean areal density for the transect and transect
        % length
        [rho_j, L_j] = calcTransectDensity(nasc(intervals,:));
        % and store
        biomass.strata(k).transect(j).name = transects(j);
        biomass.strata(k).transect(j).krillDensity = rho_j; % [g/m^2]
        biomass.strata(k).transect(j).length = L_j; % [nmi]
    end
    
    % Now calculate values per stratum
    
    % Extract the per transect data from the biomass structure for
    % convenience into variables with the same names as per EMM-16/38.
    L_j = [biomass.strata(k).transect.length]; % [nmi]
    rho_j = [biomass.strata(k).transect.krillDensity]; % [g/m^2]
    N_k = length(transects); % [1]
    
    % Compute normalized transect weighting factors
    w_j = L_j / mean(L_j); % equation 11 [1]
        
    % Compute mean krill biomass density for stratum
    %rho_k = sum(w_j .* rho_j) / sum(w_j); % as in calcDensity script
    rho_k = sum(w_j .* rho_j) / N_k; % equation 14 [g/m^2]

    % Variance component
    VarComp_k = w_j.^2 .* (rho_j - rho_k).^2; % equation 13 [g^2/m^4]

    % Compute variance of mean krill biomass density, equation 15 [g^4/m^8??]
    var_rho_k = N_k / (N_k-1) ...
        * sum(w_j.^2 .* (rho_k-rho_j).^2) / (sum(w_j)).^2 ...
        * sum(w_j.^2 .* (rho_k-rho_j).^2) / (N_k * (N_k-1));

    CV_k = 100 * sqrt(var_rho_k) / rho_k; % equation 16 [??]
    
    % Find the stratum area
    a_i = find(strcmp({strata_area.name}, strata(k)));
    
    % and store the above values in a nicer form
    biomass.strata(k).area = strata_area(a_i).area; % [km^2]
    biomass.strata(k).stratumLength = sum(L_j);
    biomass.strata(k).meanDensity = rho_k;
    biomass.strata(k).densityVariance = var_rho_k;
    biomass.strata(k).varianceComponent = VarComp_k;
    biomass.strata(k).CV = CV_k;
end

% Then calculate results for the entire survey.
% There are multiple surveys, so define which
% strata combine for each 'survey'. 

biomass.survey(1).name = 'Area 48 2019';
biomass.survey(1).strata = ["AP", "SS", "SSI", "SOI", "Sand", "ESS"];

biomass.survey(2).name = 'AMLR 2019';
biomass.survey(2).strata = ["Joinville", "Elephant", "West", "South"];

for s = 1:length(biomass.survey)
    % identify which strata are in this survey
    k = find(ismember([biomass.strata.name],biomass.survey(s).strata));
    
    % extract some per strata numbers into variables that match those used
    % in the equations
    A_k = [biomass.strata(k).area]; % [km^2]
    rho_k = [biomass.strata(k).meanDensity]; % [g/m^2]
    var_rho_k = [biomass.strata(k).densityVariance]; % [?]
    VarComp_k = [biomass.strata(k).varianceComponent]; % [?]
    
    % overall mean areal krill biomass density.
    rho = sum(A_k.*rho_k) / sum(A_k); % equation 17 [g/m^2]
    
    % overall survey variance
    var_rho = sum(A_k.^2 .* var_rho_k.^2) / sum(A_k.^2); % equation 18 [?]
    
    % overall CV of mean areal krill biomass density
    CV_rho = 100 * sqrt(var_rho) / rho; % equation 19
    
    % Biomass. 
    B_0 = sum(A_k .* rho_k); % equation 20 [km^2 * g/m^2] = [tonnes]
    
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



