%%
% Estimate the krill biomass per strata and per survey

do_define_directories

% What strata make up a survey...
surveys = get_surveys();

% Pull in the strata areas
strata_area = get_strata_areas();

% Pull in a function that converts krill length into weight
create_l_to_w_function

% Krill TS values
[krill_len, krill_sigma] = get_krill_sigma(resultsDir);

% Lengths per strata/cluster
load(fullfile(resultsDir, 'Trawls - data'), 'lf')

% Weighted conversion factor for each strata and cluster
for t = ["strata" "cluster"]
    for i = 1:length(lf.(t))
        if t == "strata"
            histedges = lf.(t)(i).ASAM2019_normalised_lf_len(1:end-1) * 1e-3; % [m]
        else
            histedges = lf.(t)(i).histedges(1:end-1) * 1e-3; % [m]
        end

        j = ismember(krill_len, histedges);
        sigma = krill_sigma(j); % [m^2]
        
        % Weighted conversion factor as per EMM-16/38, eqn 8.
        if t == "strata"
            C = sum(lf.(t)(i).ASAM2019_normalised_lf_prop .* w(histedges)) ./ ...
                sum(lf.(t)(i).ASAM2019_normalised_lf_prop .* sigma); % [g/m^2]
        else
            C = sum(lf.(t)(i).histcounts .* w(histedges)) ./ ...
                sum(lf.(t)(i).histcounts .* sigma); % [g/m^2]
        end
        
        
        lf.(t)(i).C = C;
    end
end

% for the maps
strata = jsondecode(fileread(fullfile(baseDir, repoDir, 'map_data', 'survey strata.geojson')));

% Calculate the full survey biomass. This uses the swarm detection method.
load(fullfile(resultsDir, 'NASC - data'), 'nasc')
results = runBiomass(nasc, lf, strata_area, surveys);
save(fullfile(resultsDir, 'Final results - swarm'), 'results')
writetable(results.nasc, fullfile(resultsDir, 'NASC for cross-checking.csv'))

do_areal_density_maps(results.nasc, strata, 'Krill density', resultsDir)
do_areal_density_maps(results.nasc_day, strata, 'Day krill density', resultsDir)

% Day/night differences
disp("Proportion of NASC data collected during the day: " ...
    + num2str(height(results.nasc_day)/height(results.nasc), '%.2f'))
results.biomass_strata(:,{'name' 'dayRatio'})

% Show all of the ratios, for each stratum and survey
format bank
t = [results.biomass_strata(:,{'name' 'dayRatio'}); results.biomass_survey(:,{'name' 'dayRatio'})]
format short

% Calculate a 'biomass' from the dB difference data from just KPH.
load(fullfile(resultsDir, 'NASC - data - dB difference - KPH'), 'nasc')
nasc.NASC = nasc.NASC_per_stratum;
results = runBiomass(nasc, lf, strata_area, surveys);
save(fullfile(resultsDir, 'Final results  - dB difference - KPH'), 'results')

