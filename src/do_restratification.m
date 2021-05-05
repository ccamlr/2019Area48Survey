%%
% code to take the existing echo-integrated data from the 2019 survey and 
% apply it to some different strata in several different ways. 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load in data
do_define_directories

% read in the various strata boundaries.
boundaryFiles = {'fishery_areas.geojson', 'on_off_shelf_areas.geojson', 'subareas.geojson'};

clear strata surveyed
s_i = 0;
ss_i = 0;
for i = 1:length(boundaryFiles)
    disp(['Loading ' boundaryFiles{i}])
    s = jsondecode(fileread(fullfile(baseDir, repoDir, 'map_data', boundaryFiles{i})));    

    % and pull out the useful data
    for j = 1:length(s.features)
        disp(['  processing ' s.features(j).properties.Name])
        pshape = polyshape();
        if iscell(s.features(j).geometry.coordinates)
            for k = 1:length(s.features(j).geometry.coordinates)
                if ~iscell(s.features(j).geometry.coordinates{k})
                    pshape = addboundary(pshape, squeeze(s.features(j).geometry.coordinates{k}));
                else
                    for kk = 1:length(s.features(j).geometry.coordinates{k})
                        pshape = addboundary(pshape, s.features(j).geometry.coordinates{k}{kk});
                    end
                end
            end
        else
           pshape = polyshape(squeeze(s.features(j).geometry.coordinates)); 
        end
        temp_strata.name = s.features(j).properties.Name;
        temp_strata.shape = pshape;
        temp_strata.area = s.features(j).properties.area;
        
        if ~endsWith(s.features(j).properties.Name, 'surveyed')
            s_i = s_i + 1;
            strata(s_i) = temp_strata;
        else
            ss_i = ss_i + 1;
            surveyed(ss_i) = temp_strata;
        end
        
    end
end

% read in the krill areal density values
load(fullfile(resultsDir, 'Final results - swarm'), 'results')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% make a map of what we have
clf
m_proj('Azimuthal Equal-area', 'long', -41.5, 'lat', -61.2, ...
    'radius', 20, 'rot', 0, 'rectbox', 'on');
m_gshhs_i('patch', 'FaceColor', [0.5 0.5 0.5]); % the land
hold on

% the existing backscatter data
m_plot(results.nasc.Longitude, results.nasc.Latitude, '.')

% the strata boundaries
for i = 1:length(strata)
    [x,y] = boundary(strata(i).shape);
    m_line(x, y, 'color', 'r')
end

m_grid()

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculate biomasses

nasc_all = results.nasc;
% remove nasc points where the krill density is NaN. These are areas were
% we didn't have a krill distribution to calculate the conversion factor.
i = isnan(nasc_all.rho);
%unique(nasc_all(i,:).Stratum); % seems to just be SA483 stratum
nasc_all = nasc_all(~i,:);

% remove nasc points from the AMLR transects
AMLR = {'Bransfield', 'Elephant', 'Joinville', 'West'};
to_keep = true(height(nasc_all),1);
for i = 1:length(AMLR)
    to_keep = to_keep & ~strcmp(nasc_all.Stratum, AMLR{i});
end
nasc_all = nasc_all(to_keep,:);

earth_radius = 6370e3; % [m]
 
biomass_results = struct([]);

for i = 1:length(strata)
    disp(['Processing ' strata(i).name])
    
    % find all nasc points that are within the current stratum. 
    s = isinterior(strata(i).shape, nasc_all.Longitude, nasc_all.Latitude);
    % And pull out into a variable for convenience
    nasc = nasc_all(s,:); 
    
    transects = unique(nasc.Transect);
    
    % get the area of the surveyed part of the stataum. These have the same
    % name, but have 'surveyed' appended.
    f = find(strcmp({surveyed.name}, [strata(i).name ' surveyed']));    
    surveyed_area = surveyed(f).area * 1e6; % [m^2]
    stratum_area = strata(i).area * 1e6; % [m^2]
 
    % now do a mean krill density in each transect in the current stratum
    % and calculate the transect length.
    transect_rho = 0.0;
    transect_lengths = [];
    for t_i = 1:length(transects)
        k = strcmp(nasc.Transect, transects(t_i));
        transect_rho(t_i) = mean(nasc.rho(k) * 1e-3); % [kg/m^2]
        transect_lengths(t_i) = sum(k) * 1852; % [m] assuming each NASC value represents 1 n.mile
    end
    
    % and then the mean rho for the stratum, weighted by transect length
    w = transect_lengths ./ (1.0/length(transect_lengths) * sum(transect_lengths));
    surveyed_rho = 1.0 / length(w) * nansum(w.* transect_rho); % [kg/m^2]
    
    surveyed_biomass = surveyed_rho * surveyed_area * 1e-3; % [tonnes]
    stratum_biomass  = surveyed_rho * stratum_area * 1e-3; % [tonnes]
 
    biomass_results(i).stratum_name = strata(i).name;
    biomass_results(i).surveyed_rho = surveyed_rho * 1000; % [g/m^2]
    biomass_results(i).surveyed_area_biomass = surveyed_biomass; % [t]
    biomass_results(i).stratum_biomass = stratum_biomass; % [t]
    biomass_results(i).surveyed_area = surveyed_area * 1e-6; % [km^2]
    biomass_results(i).stratum_area = stratum_area * 1e-6; % [km^2]
    biomass_results(i).area_ratio = stratum_area / surveyed_area; % [1]
    
    % do a map for each strata
    clf
    m_proj('Azimuthal Equal-area', 'long', -41.5, 'lat', -61.2, ...
           'radius', 20, 'rot', 0, 'rectbox', 'on');
    m_gshhs_i('patch', 'FaceColor', [0.5 0.5 0.5]); % the land

    % the strata boundary
    [x,y] = boundary(strata(i).shape);
    m_line(x, y, 'color', 'r')
    hold on
    % the surveyed boundary
    [x,y] = boundary(surveyed(f).shape);
    m_line(x, y, 'color', 'g')

    % the backscatter data that is being used
    m_plot(nasc.Longitude, nasc.Latitude, '.')
 
    title([strata(i).name])
    m_grid()

    % save the map
    print(fullfile(resultsDir, ['Restratification ' strata(i).name '.png']), '-r300', '-dpng');
end

biomass_results_table = struct2table(biomass_results);
writetable(biomass_results_table, fullfile(resultsDir, 'Restratification biomass.csv'))

% There are 3 ways we'll re-stratify:
% a) Use the transects within a given subarea to estimate the biomass with
% the area surveyed in that subarea. Then extrapolate to the area of the
% subarea using the survey mean density. 
%
% b) Use the transects within a given subarea to estimate the biomass
% within the area surveyed, but subdivide each survey into on-shelf and
% off-shelf strata. Then extrapolate to the area of on-shelf and
% off-shelf in the subarea using the survey strata means (g/m2).  
%
% c) Use the transects within a given subarea to estimate the biomass
% within the area usually fished in that subarea, based only on the survey
% transects in that fished area and extrapolated to the fished area using
% the mean (g/m2).

% 
