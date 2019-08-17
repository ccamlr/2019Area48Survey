function [results, nasc, nasc_day] = runBiomass(nasc, lf, strata_area, surveys)

% Runs the biomass calculations from the given nasc values, trawl lf's,
% strata areas and surveys

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

% Where the real numbers get made...
results = calcBiomass(nasc, strata_area, surveys);

% and what if we exclude the night data?
nasc_day = nasc(nasc.civilDaytime == true,:);
results_day = calcBiomass(nasc_day, strata_area, surveys);

% and put the daytime density ratios (and biomass, it's the same number) into the results structure
results.biomass_strata.dayRatio = results_day.biomass_strata.meanDensity ./ results.biomass_strata.meanDensity;
results.biomass_survey.dayRatio = results_day.biomass_survey.meanDensity ./ results.biomass_survey.meanDensity;
results.nasc = nasc;
results.nasc_day = nasc_day;


