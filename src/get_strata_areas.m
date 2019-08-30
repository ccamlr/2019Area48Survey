function strata_area = get_strata_areas()

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
    'Joinville', 18151; 'SOF', 0; 'SOC', 0; 'WCB', 0};

% and turn the areas into a more useful structure
strata_area(size(s,1)) = struct('name',[], 'area', []);
for i = 1:size(s,1)
    strata_area(i).name = s{i,1};
    strata_area(i).area = s{i,2};
end

