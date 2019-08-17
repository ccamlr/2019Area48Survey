function results = calcBiomass(nasc, strata, surveys)

    % Now calculate biomass density for each transect, stratum and finally
    % survey. This follows the equations in EMM-16/38 section 9, with some
    % corrections.
    
    % We are given a list of strata to work with, but sometimes we have no
    % data for some of the strata, so deal with that where necessary.

    for k = 1:length(strata) % loop over all strata
        results.strata(k).name = strata(k).name;
        transects = unique(nasc.Transect(nasc.Stratum == strata(k).name));
        for j = 1:length(transects) % loop over all transects in each strata
            % Find all NASC values for the current transect
            intervals = find(nasc.Transect == transects(j) & nasc.Stratum == strata(k).name);
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
        if isempty(results.strata(k).transect) % no data for this stratum
            L_j = 0;
            rho_j = NaN;
            N_k = 0;
        else
            L_j = [results.strata(k).transect.length]; % [nmi]
            rho_j = [results.strata(k).transect.krillDensity]; % [g/m^2]
            N_k = length(transects); % [1]
        end
      
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
        
        % and store the above values in a nicer form
        results.strata(k).area = strata(k).area; % [km^2]
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
    for s = 1:length(surveys)
        results.survey(s).name = surveys(s).name;
        results.survey(s).strata = surveys(s).strata;
        
        % identify which strata are in this survey
        k = find(ismember({results.strata.name}, surveys(s).strata));
        
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
    t.Properties.VariableNames{end} = 'varianceComponent';
    t.Properties.VariableUnits={'' 'km^2' 'gm^-2' 't' 't^2'};
    results.biomass_strata = t;
end

