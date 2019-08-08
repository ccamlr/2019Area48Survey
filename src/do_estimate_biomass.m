%%
% Estimate the krill biomass per strata and per survey

baseDir = 'I:\KRILL2019';
repoDir = '2019Area48SurveyRepo';
resultsDir = fullfile(baseDir, repoDir, 'results');

% Strata areas. 

% ** Should be calculated directly from the latest strata boundaries **

% CCAMLR 2000 strata areas from Bo workshop 2004, section 2.3
% AMLR strata areas from Figure 1 of EMM-11/26 and confirmed via email from
% Christian Reiss. He notes that for Joinville they tend not to report
% biomass (only density) as their survey coverage is quite variable from year to year. 

% Units are km^2
s = {'AP', 473318; 'SS' 1109789; 'ESS', 321800; ...
    'SSI', 48654; 'SOI', 24409; 'SG', 25000; 'Sand', 62274; ...
    'Joinville', 18151; 'Elephant', 43865; 'West', 38524; 'Bransfield', 24479; ...
    'SOF', 0; 'SOC', 0};
% and turn this into a more useful structure
strata_area(size(s,1)) = struct('name',[], 'area', []);
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

% and to produce plots of krill areal density later on, calculate a krill
% density for each nasc value.
nasc.rho = nasc.NASC .* nasc.C / 1852^2;

% and set the units for the nasc columns to make it clear what is what.
% Uses the UDUNITS convention.
nasc.Properties.VariableUnits = {'' 'degrees_north' 'degrees_east' ...
    'm^2 nautical_mile^-2' '' '' '' 'g m^-2' 'g m^-2'};

% Now calculate biomass density for each transect, stratum and finally
% survey. This follows the equations in EMM-16/38 section 9, with some
% corrections.

clear results
for k = 1:length(strata) % loop over all strata 
    results.strata(k).name = strata(k);
    transects = unique(nasc.Transect(nasc.Stratum == strata(k)));
    for j = 1:length(transects) % loop over all transects in each strata
        % Find all NASC values for the current transect
        intervals = find(nasc.Transect == transects(j) & nasc.Stratum == strata(k));
        % Calculate the mean areal density for the transect and transect
        % length
        [rho_j, L_j] = calcTransectDensity(nasc(intervals,:));
        % and store
        results.strata(k).transect(j).name = transects(j);
        results.strata(k).transect(j).krillDensity = rho_j; % [g/m^2]
        results.strata(k).transect(j).length = L_j; % [nmi]
    end
    
    % Now calculate values per stratum
    
    % Extract the per transect data from the results structure for
    % convenience into variables with the same names as per EMM-16/38.
    L_j = [results.strata(k).transect.length]; % [nmi]
    rho_j = [results.strata(k).transect.krillDensity]; % [g/m^2]
    N_k = length(transects); % [1]
    
    % Compute normalized transect weighting factors
    w_j = L_j / mean(L_j); % equation 11 [1]
        
    % Compute mean krill biomass density for stratum
    %rho_k = sum(w_j .* rho_j) / sum(w_j); % as in calcDensity script
    rho_k = sum(w_j .* rho_j) / N_k; % equation 14 [g/m^2]

    % Variance component
    VarComp_j = w_j.^2 .* (rho_j - rho_k).^2; % equation 13 [g^2/m^4]

    % Variance of mean krill biomass density, modified version of equation
    % 15 that follows the 'where' part of eqn 3 in Jolly & Hampton 1990
    % [g^2/m^4]
    var_rho_k = sum(w_j.^2 .* (rho_k-rho_j).^2) / (N_k * (N_k-1));

    CV_k = 100 * sqrt(var_rho_k) / rho_k; % equation 16 [1]
    
    % Find the stratum area
    a_i = find(strcmp({strata_area.name}, strata(k)));
        
    % and store the above values in a nicer form
    results.strata(k).area = strata_area(a_i).area; % [km^2]
    results.strata(k).stratumLength = sum(L_j);
    results.strata(k).meanDensity = rho_k;
    results.strata(k).densityVariance = var_rho_k;
    results.strata(k).varianceComponent_transect = VarComp_j;
    results.strata(k).CV = CV_k;
    results.strata(k).biomass = results.strata(k).area .* results.strata(k).meanDensity; % [t]
    
    % variance component for the stratum. Not given in EMM-16/38.
    VarComp_k = results.strata(k).area.^2 .* var_rho_k; % [t^2]
    results.strata(k).varianceComponent_stratum = VarComp_k;

end

% Then calculate results for the entire survey.
% There are multiple surveys, so define which
% strata combine for each 'survey'. 

results.survey(1).name = 'Area 48 2019';
results.survey(1).strata = ["AP", "SS", "SSI", "SOI", "Sand", "ESS"];

results.survey(2).name = 'AMLR 2019';
results.survey(2).strata = ["Joinville", "Elephant", "West", "Bransfield"];

results.survey(3).name = 'Norway1';
results.survey(3).strata = ["SOF"];

results.survey(4).name = 'Norway1';
results.survey(4).strata = ["SOC"];


for s = 1:length(results.survey)
    % identify which strata are in this survey
    k = find(ismember([results.strata.name],results.survey(s).strata));
    
    % extract some per strata numbers into variables that match those used
    % in the equations
    A_k = [results.strata(k).area]; % [km^2]
    rho_k = [results.strata(k).meanDensity]; % [g/m^2]
    var_rho_k = [results.strata(k).densityVariance]; % [g^2/m^4]
    VarComp_k = [results.strata(k).varianceComponent_stratum]; % [t]
            
    % overall mean areal krill biomass density.
    rho = sum(A_k.*rho_k) / sum(A_k); % equation 17 [g/m^2]
    
    % overall survey variance. Modified version of equation 18 that now
    % follows eqn 3 of Jolly & Hampton 1990.
    var_rho = sum(A_k.^2 .* var_rho_k) / sum(A_k)^2; % [g^2/m^4]
    
    % overall CV of mean areal krill biomass density
    CV_rho = 100 * sqrt(var_rho) / rho; % equation 19 [%]
    
    % Biomass. 
    B_0 = sum(A_k .* rho_k); % equation 20 [km^2 * g/m^2] = [tonnes]
    
    % Variance of B_0
    var_B_0 = sum(VarComp_k); % equation 21 [t^2]
    
    % Overall coefficient of variation of B_0. Equation 22 [1].
    CV_B_0 = 100 * sqrt(var_B_0) ./ B_0; % [%]
    
    % and store the results in a nicer form
    results.survey(s).meanDensity = rho; % [g m^-2]
    results.survey(s).meanDensityVariance = var_rho; % [g^2 m^-4]
    results.survey(s).meanDensityCV = CV_rho; % [%]

    results.survey(s).biomass = B_0; % [t]
    results.survey(s).variance = var_B_0; % [t^2]
    results.survey(s).CV = CV_B_0; % [%]
end

% now put some per strata and per survey data into a table for nicely
% formatted output

t = struct2table(results.survey);
t = removevars(t, {'strata'});
t.Properties.VariableUnits = {'','gm^-2','g^2m^-4','%','t','t^2','%'};
results.biomass_survey = t;

% Follow the contents and arrangement of Table 4 in EMM-11/20
t = struct2table(results.strata);
t = removevars(t, {'transect', 'stratumLength' 'CV' 'densityVariance' 'varianceComponent_transect'});
t = movevars(t,'biomass','Before','varianceComponent_stratum');
t = t([2 1 7 3 4 6 5 9 10 11 8 12 13],:); % adhocery... sort rows to match some existing reports
t.Properties.VariableNames{end} = 'varianceComponent';
t.Properties.VariableUnits={'' 'km^2' 'gm^-2' 't' 't^2'};
results.biomass_strata = t;

% and save the results
save(fullfile(resultsDir, 'Final results'), 'results')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% And now produce some plots of krill areal density. These are very similar
% to what is done in the do_combine_nasc.m script, but use areal density
% rather than NASC.

strata = jsondecode(fileread(fullfile(baseDir, repoDir, 'map_data', 'survey strata.geojson')));

% Use the same max scale across all plots
maxRho = max(nasc.rho);
maxSize = 200; % [points^2] of drawn circles
legendScatterSizes = [50 500 2500 5000]; % [g/m^2]

% Map coloured by vessel
figure(1)
clf
plot_standard_map(strata, 'showStrataNames', false)

% Use a different colour for each vessel
v = unique(nasc.Vessel);
h = nan(length(v), 1);
for i = 1:length(v)
    j = find(nasc.Vessel == v(i));
    h(i) = m_scatter(nasc.Longitude(j), nasc.Latitude(j), nasc.rho(j)/maxRho*maxSize+1, 'filled');
end

m_grid('box', 'on')
legend(h, v, 'Location', 'SouthEast')

print(fullfile(resultsDir, 'Krill density - by vessel'), '-dpng','-r300')

% Map coloured by stratum. Too crowded to be really useful...
figure(2)
clf
plot_standard_map(strata, 'showStrataNames', false)

% Use a different colour and symbol for each statum

s = unique(nasc.Stratum);
symbols = {'o' 'o' 'o' 'o' 'o' 'o' 'o' 'd' 'd' 'd' 'd' 'd' 'd' 'd'};
h = [];
for i = 1:length(s)
    j = find(nasc.Stratum == s(i));
    h(i) = m_scatter(nasc.Longitude(j), nasc.Latitude(j), nasc.rho(j)/maxRho*maxSize+1, 'filled', symbols{i});
end

m_grid('box', 'on')
legend(h, s, 'Location', 'SouthEast', 'NumColumns', 2, 'Interpreter', 'none')

print(fullfile(resultsDir, 'Krill density - by stratum'), '-dpng','-r300')

%%%%%%%%%%%
figure(3)
clf

s = ["Bransfield" "Elephant" "Joinville" "West"];
plot_standard_map(strata, 'centrePoint', [-58 -62], 'radius', 4, ...
    'strata', s, 'showStrataNames', false, ...
    'coastDetail', 'high')

for i = 1:length(s)
    j = find(nasc.Stratum == s(i));
    m_scatter(nasc.Longitude(j), nasc.Latitude(j), nasc.rho(j)/maxRho*maxSize+1, 'filled', 'o');
end
plot_standard_map_rho_legend(legendScatterSizes, maxRho, maxSize)

print(fullfile(resultsDir, 'Krill density - AMLR'), '-dpng','-r300')

%%%%%%%%%%%
figure(4)
clf

s = ["ESS" "Sand" "SG" "SS" "AP" "SSI" "SOI"];
plot_standard_map(strata, 'centrePoint', [-45 -60], 'radius', 17.5, ...
    'strata', s, 'showStrataNames', false, ...
    'coastDetail', 'intermediate')

for i = 1:length(s)
    j = find(nasc.Stratum == s(i));
     m_scatter(nasc.Longitude(j), nasc.Latitude(j), nasc.rho(j)/maxRho*maxSize+1, 'filled', 'o');
end
plot_standard_map_rho_legend(legendScatterSizes, maxRho, maxSize)

print(fullfile(resultsDir, 'Krill density - CCAMLR 2000'), '-dpng','-r300')

%%%%%%%%%%%
figure(5)
clf

s = ["SOI" "SOC" "SOF"];
plot_standard_map(strata, 'centrePoint', [-45.7 -60.75], 'radius', 2.5, ...
    'strata', s, 'showStrataNames', false, ...
    'coastDetail', 'fine')

for i = 1:length(s)
    j = find(nasc.Stratum == s(i));
    m_scatter(nasc.Longitude(j), nasc.Latitude(j), nasc.rho(j)/maxRho*maxSize+1, 'filled', 'o');
end
plot_standard_map_rho_legend(legendScatterSizes, maxRho, maxSize)

print(fullfile(resultsDir, 'Krill density - South Orkney'), '-dpng','-r300')




