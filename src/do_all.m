%%
% This script lists all the scripts neeed to calculate krill biomass, in a
% sensible order.

% Converts ek80 .raw files into smaller files and converts LSSS work files
% into Echoview line and region files.
do_make_ready_for_echoview

% The next step here is to run the EchoviewR_integrate.R script to
% echo-integrate the acoustic data. Two vessels are not run by that script,
% as the individual countries do that themselves (UK and China).

% Work out sound speed and absorptions to use
do_average_ctds

% Work out krill lfs to use
do_cluster_trawls

% Combine nasc outputs from each ship
do_combine_nasc

% Calculate krill TS using the survey parameters.
% This only needs to be done once - the output is left in the 'results'
% directory and is used by other code as necessary.
do_calculate_ts

%% Calculate the biomass
do_estimate_biomass

%% Work up the data we have where the 3-freq dB-difference method is pissible
do_db_difference

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Show the biomass results 
do_define_directories
load(fullfile(resultsDir, 'Final results - swarm'), 'results')

% Show the per strata results
disp('Results per strata:')
disp(results.biomass_strata)

% Show the per survey results
disp('Results per survey:')
disp(results.biomass_survey)

% Do it to make it easy to copy/paste an appropriately formatted
% table into the ASAM report (via Excel)
disp('strata for pasting:')

jf = java.text.DecimalFormat; % formatting with thousands commas
for i = 1:height(results.biomass_strata)
    r = results.biomass_strata(i,{'area' 'meanDensity' 'biomass' 'varianceComponent'});
    fprintf('%s; %.1f; %s; %s,000\n', jf.format(r.area), r.meanDensity, ...
        jf.format(round(r.biomass/1e3)*1e3), ...
        jf.format(round(r.varianceComponent/1e9)))
end

% surveys
disp('surveys for pasting')
for i = 1:height(results.biomass_survey)
    r = results.biomass_survey(i,{'meanDensity' 'meanDensityVariance' 'meanDensityCV' 'biomass' 'variance' 'CV'});
    fprintf('%.1f; %.1f; %.0f; %s; %s,000; %.0f\n', r.meanDensity, ...
        r.meanDensityVariance, r.meanDensityCV, ...
        jf.format(round(r.biomass/1e3)*1e3), ...
        jf.format(round(r.variance/1e9)), ...
        r.CV)
end

% data only biomass per strata
disp('day only strata biomass')
for i = 1:height(results.biomass_strata_day)
    r = results.biomass_strata_day(i,{'name' 'meanDensity' 'biomass'});
    fprintf('%s; %.1f; %s\n', r.name{1}, ...
        r.meanDensity, ...
        jf.format(round(r.biomass/1e3)*1e3))
end

% report the start and stop survey dates for each stratum
strata = unique(results.nasc.Stratum);
for i = 1:length(strata)
    s = strcmp(results.nasc.Stratum, strata(i));
    startTime = min(results.nasc.Ping_timestamp(s));
    stopTime = max(results.nasc.Ping_timestamp(s));
    fprintf('%s, %s, %s, %d\n', strata(i), datestr(startTime, 'dd mmm yyyy'), ...
        datestr(stopTime, 'dd mmm yyyy'), sum(s))
end


%% Output the results in the CCAMLR-formatted spreadsheet format
% This is format for WG-ASAM to use to collate acoustic survey results
% Use the generated csv file to import into the spreadsheet

t = struct2table(results.strata);
t = removevars(t, {'transect', 'stratumLength', 'densityVariance', 'varianceComponent_transect', 'biomass', 'varianceComponent_stratum'});
t.Properties.VariableNames{'name'} = 'stratumName';
t.Properties.VariableNames{'CV'} = 'densityCV';
t.Properties.VariableNames{'area'} = 'surveyArea';
t.Properties.VariableNames{'meanDensity'} = 'densityEstimate';

t.year = repmat(2019, height(t), 1);
t.month = repmat({'Dec,Jan,Feb,Mar'}, height(t), 1);
t.vessel = repmat({'several'}, height(t), 1);
t.contributor = repmat({'ARK,China,Korea,Norway,Ukraine,United Kingdom'}, height(t), 1);
t.subarea = repmat({''}, height(t), 1);
t.surveyName = repmat({'2019 synoptic'}, height(t), 1);
t.CVmethod = repmat({'Jolly and Hampton'}, height(t), 1);
t.echosounder = repmat({'EK60,ES70,ES80,EK80'}, height(t), 1);
t.frequency = repmat('120', height(t), 1);
t.otherFreq = repmat({'various'}, height(t), 1);
t.targetIDMethod = repmat({'Swarms identification'}, height(t), 1);
t.dbDifference = repmat({'NA'}, height(t), 1);
t.TSmodel = repmat({'full SDWBA (ASAM 2010 parameterisation)'}, height(t), 1);
t.depthRange = repmat(250, height(t), 1);
t.timeOfDay = repmat({'day and night'}, height(t), 1);
t.surveyDesign = repmat({'As per CCAMLR 2000 survey'}, height(t), 1);
t.reference = repmat({''}, height(t), 1);
t.note = repmat({''}, height(t), 1);

t = movevars(t, {'densityEstimate'}, 'After', 'surveyName');
t = movevars(t, {'densityCV'}, 'After', 'densityEstimate');
t = movevars(t, {'surveyArea'}, 'After', 'CVmethod');
t = movevars(t, {'stratumName'}, 'After', 'timeOfDay');

% fix some things
for i = 1:height(t)
    switch t.stratumName{i}
        case {'SS'}
            t.subarea(i) = {'48.2, 48.3'};
        case {'AP', 'SSI', 'Elephant', 'West', 'Bransfield', 'Joinville'}
            t.subarea(i) = {'48.1'};
        case {'SOI', 'SOF', 'SOC'}
            t.subarea(i) = {'48.2'};
        case {'SG', 'WCB'}
            t.subarea(i) = {'48.3'};
        case {'Sand', 'ESS'}
            t.subarea(i) = {'48.4'};
    end
    
    %if strcmp(t.stratumName{i}, 'WCB')
    %    t.contributor(i) = {'United Kingdom'};
    %end
end



writetable(t, fullfile(resultsDir, 'WG_EMM_krill_biomass_time_series_information_synoptic_2019.csv'), ...
    'FileType', 'text', 'Delimiter', ',', 'QuoteStrings', true)

%% Show the comparison of KPH swarm verses dB-difference krill classification
do_define_directories
resultsDir = resultsDir2019; % Use old data for this.
load(fullfile(resultsDir, 'Final results - swarm'), 'results')
results_swarm = results;
load(fullfile(resultsDir, 'Final results  - dB difference - KPH'), 'results')
results_dBdiff = results;

% Mainly interested in the density per transect, and per-strata biomass, so
% do comparison tables of that.
clear transect
k = 1;
for i = 1:length(results_dBdiff.strata)
    if ~isempty(results_dBdiff.strata(i).transect)
        for j = 1:length(results_dBdiff.strata(i).transect)
            name = string(results_dBdiff.strata(i).name) + results_dBdiff.strata(i).transect(j).name;
            transect.name(k,1) = name;
            transect.strata(k,1) = string(results_dBdiff.strata(i).name);
            transect.rho_dB(k,1) =  results_dBdiff.strata(i).transect(j).krillDensity;
            transect.rho_swarm(k,1) = results_swarm.strata(i).transect(j).krillDensity;
            transect.length(k,1) = results_swarm.strata(i).transect(j).length;
            k = k + 1;
        end
    end
end
transect = struct2table(transect);
transect.prop = transect.rho_dB ./ transect.rho_swarm;

% Per-stratum comparison table
clear strata
k = 1;
for i = 1:length(results_dBdiff.strata)
    if ~isnan(results_dBdiff.strata(i).meanDensity)
        strata.name(k,1) = string(results_dBdiff.strata(i).name);
        strata.rho_dB(k,1) = results_dBdiff.strata(i).meanDensity;
        strata.rho_swarm(k,1) = results_swarm.strata(i).meanDensity;
        strata.area(k,1) = results_dBdiff.strata(i).area;
        k = k + 1;
    end
end
strata = struct2table(strata);
strata.prop = strata.rho_dB ./ strata.rho_swarm;


% for copy/pasting into Word, etc
disp('Per-transect densities')
for i = 1:height(transect)
    disp(transect.name(i) + ' ' + num2str(transect.rho_swarm(i), '%.1f') + ' ' ...
        + num2str(transect.prop(i), '%.2f') + ' ' ...
        + num2str(transect.length(i), '%.f'))
end

disp('Per-stratum biomass comparison')
for i = 1:height(strata)
    disp(strata.name(i) + ' ' + num2str(strata.rho_swarm(i), '%.1f') + ' ' ...
        + num2str(strata.prop(i), '%.1f') + ' ' ...
        + num2str(strata.area(i), '%.f'))    
end

% A plot of the data with regression
mdl = fitlm(transect.rho_swarm, transect.rho_dB);
ci = coefCI(mdl, 0.05);
% avoids Stats toolbox
% p = polyfit(transect.rho_swarm, transect.rho_dB, 1) 
eqn = {"\rho_{dB} = (" + num2str(mdl.Coefficients.Estimate(2), '%.1f') + ...
       " \pm " + num2str(mdl.Coefficients.Estimate(2) - ci(2,1), '%.1f') + ...
       ")\rho_{swarm} + (" + num2str(mdl.Coefficients.Estimate(1), '%.1f') + ...
       " \pm " + num2str(mdl.Coefficients.Estimate(1) - ci(1,1), '%.1f)')
    "r^2 = " + num2str(mdl.Rsquared.Adjusted, '%.3f')};
m = max([transect.rho_swarm; transect.rho_dB]);

figure(1)
clf
strata_symbols = {'p', '<', 'o', 's', 'd'};
strata = unique(transect.strata);
for i = 1:length(strata)
    j = find(transect.strata == strata(i));
    plot(transect.rho_swarm(j), transect.rho_dB(j), 'LineStyle', 'none', ...
        'Marker', strata_symbols{i}, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k')
    hold on
end
plot([0 m],  mdl.Coefficients.Estimate(2) * [0 m] + mdl.Coefficients.Estimate(1), 'k', 'LineWidth', 1.5)
grid
xlabel('Areal density, swarms (g m^{-2})')
ylabel('Areal density, dB-difference (g m^{-2})')

set(gca, 'XLim', [0 m], 'YLim', [0 m])
textLoc(eqn, 'NorthWest');

ifile = fullfile(resultsDir, 'Swarms verses dB-diff - regression.png');
print(ifile, '-r300', '-dpng')
crop_image(ifile)

