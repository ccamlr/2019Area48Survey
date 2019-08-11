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
    'Joinville', 18151, 'SOF', 0; 'SOC', 0};

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
    'm^2 nautical_mile^-2' '' '' '' 'g m^-2' 'g m^-2'};

% Where the real numbers get made...
results = calcBiomass(nasc, strata_area, surveys);

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
h = nan(size(s));
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
    'strata', s, 'showStrataNames', true, ...
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
    'strata', s, 'showStrataNames', true, ...
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
    'strata', s, 'showStrataNames', true, ...
    'coastDetail', 'fine')

for i = 1:length(s)
    j = find(nasc.Stratum == s(i));
    m_scatter(nasc.Longitude(j), nasc.Latitude(j), nasc.rho(j)/maxRho*maxSize+1, 'filled', 'o');
end
plot_standard_map_rho_legend(legendScatterSizes, maxRho, maxSize)

print(fullfile(resultsDir, 'Krill density - South Orkney'), '-dpng','-r300')




