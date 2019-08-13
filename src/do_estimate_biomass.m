%%
% Estimate the krill biomass per strata and per survey

do_define_directories

% What strata make up a survey...
s = 1;
surveys(s).name = 'CCAMLR 2019';
surveys(s).strata = ["AP", "SS", "SSI", "SOI", "Sand", "ESS"]; s = s + 1;

surveys(s).name = 'AMLR 2019';
surveys(s).strata = ["Joinville", "Elephant", "West", "Bransfield"]; s = s + 1;

surveys(s).name = 'Norway1';
surveys(s).strata = "SOF"; s = s + 1;

surveys(s).name = 'Norway1';
surveys(s).strata = "SOC"; s = s + 1;

% Strata areas. 
% ** Should be calculated directly from the latest strata boundaries **

% CCAMLR 2000 strata areas from Bo workshop 2004, section 2.3
% AMLR strata areas from Figure 1 of EMM-11/26 and confirmed via email from
% Christian Reiss. He notes that for Joinville they tend not to report
% biomass (only density) as their survey coverage is quite variable from year to year. 
% Strata names and area [km^2]
s = {'AP', 473318; 'SS' 1109789; 'ESS', 321800; ...
    'SSI', 48654; 'SOI', 24409; 'SG', 25000; 'Sand', 62274; ...
    'Elephant', 43865; 'West', 38524; 'Bransfield', 24479; ...
    'Joinville', 18151; 'SOF', 0; 'SOC', 0};

% and turn the areas into a more useful structure
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

% Associate conversion factor with each NASC value

% This is the slow way, but will catch nasc rows with unknown strata...
for i = 1:height(nasc)
    j = find(strcmp(nasc.Stratum(i), {lf.strata.stratum}));
    if ~isempty(j)
        nasc.C(i) = lf.strata(j).C;
    else
        nasc.C(i) = NaN; % for nasc values in strata without a krill lf, or where a strata is not present.
    end
end

% try assigning all NASC to a single cluster lf
%nasc.C = repelem(lf.cluster(1).C, height(nasc), 1);

% A krill density for each nasc value.
nasc.rho = nasc.NASC .* nasc.C / 1852^2;
% and set the units for the nasc columns to make it clear what is what.
% Uses the UDUNITS convention.
nasc.Properties.VariableUnits = {'' 'degrees_north' 'degrees_east' ...
    'm^2 nautical_mile^-2' '' '' '' 'g m^-2' 'g m^-2', ''};

% Where the real numbers get made...
results = calcBiomass(nasc, strata_area, surveys);

% and what if we exclude the night data?
nasc_day = nasc(nasc.civilDaytime == true,:);
results_day = calcBiomass(nasc_day, strata_area, surveys);

% and put the daytime density ratios (and biomass, it's the same number) into the results structure
results.biomass_strata.dayRatio = results_day.biomass_strata.meanDensity ./ results.biomass_strata.meanDensity;
results.biomass_survey.dayRatio = results_day.biomass_survey.meanDensity ./ results.biomass_survey.meanDensity;

% and save the results
save(fullfile(resultsDir, 'Final results'), 'results')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% And now produce some plots of krill areal density. These are very similar
% to what is done in the do_combine_nasc.m script, but use areal density
% rather than NASC.
strata = jsondecode(fileread(fullfile(baseDir, repoDir, 'map_data', 'survey strata.geojson')));
do_areal_density_maps(nasc, strata, 'Krill density', resultsDir)
do_areal_density_maps(nasc_day, strata, 'Day krill density', resultsDir)


