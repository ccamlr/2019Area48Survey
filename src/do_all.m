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

%% Show the biomass results 
do_define_directories
load(fullfile(resultsDir, 'Final results'), 'results')

% Show the per strata results
disp('Results per strata:')
disp(results.biomass_strata)

% Show the per survey results
disp('Results per survey:')
disp(results.biomass_survey)


% Do it to make it easy to copy/paste an appropriately formatted
% table into the ASAM report (via Excel)
disp('strata for pasting:')

for i = 1:height(results.biomass_strata)
    r = results.biomass_strata(i,2:5);
    fprintf('%s; %.1f; %s; %s,000\n', jf.format(r.area), r.meanDensity, ...
        jf.format(round(r.biomass/1e3)*1e3), ...
        jf.format(round(r.varianceComponent/1e9)))
end

% surveys
disp('surveys for pasting')
for i = 1:height(results.biomass_survey)
    r = results.biomass_survey(i,2:7);
    fprintf('%.1f; %.1f; %.0f; %s; %s,000; %.0f\n', r.meanDensity, ...
        r.meanDensityVariance, r.meanDensityCV, ...
        jf.format(round(r.biomass/1e3)*1e3), ...
        jf.format(round(r.variance/1e9)), ...
        r.CV)
end
