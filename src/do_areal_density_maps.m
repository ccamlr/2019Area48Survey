function do_areal_density_maps(nasc, strata, prefix, saveDir)
% Make various plots of the krill areal density. 

    % Use the same max scale across all plots
    maxRho = max(nasc.rho);
    maxSize = 200; % [points^2] of drawn circles
    legendScatterSizes = [50 500 5000]; % [g/m^2]

    % Map coloured by vessel, CCAMLR 2000 strata
    figure(1)
    clf
    s = ["ESS" "Sand" "SG" "SS" "AP" "SSI" "SOI"];
    plot_standard_map(strata, 'showStrataNames', false, 'strataToPlot', s, 'rotation', 0)

    % Use a different colour for each vessel
    %v = unique(nasc.Vessel);
    v = ["KPH" "CDH" "MS" "RRS"];
    h = nan(length(v), 1);
    for i = 1:length(v)
        j = find(nasc.Vessel == v(i));
        h(i) = m_scatter(nasc.Longitude(j), nasc.Latitude(j), nasc.rho(j)/maxRho*maxSize+1, 'filled');
    end

    m_grid('box', 'on')
    legend(h, v, 'Location', 'SouthEast')

    ifile = fullfile(saveDir, [prefix ' - by vessel - CCAMLR 2000 strata.png']);
    print(ifile, '-dpng','-r300')
    crop_image(ifile)
    
    % Map coloured by vessel, AMLR strata
    figure(2)
    clf
    s = ["Bransfield" "Elephant" "Joinville" "West"];

    plot_standard_map(strata, 'centrePoint', [-58 -62], 'radius', 4, ...
        'strataToPlot', s, 'showStrataNames', true, ...
        'coastDetail', 'high')

    % Use a different colour for each vessel
    %v = unique(nasc.Vessel);
    v = ["FRH" "KJH"];
    h = nan(length(v), 1);
    for i = 1:length(v)
        j = find(nasc.Vessel == v(i));
        h(i) = m_scatter(nasc.Longitude(j), nasc.Latitude(j), nasc.rho(j)/maxRho*maxSize+1, 'filled');
    end

    m_grid('box', 'on')
    legend(h, v, 'Location', 'SouthEast')

    ifile = fullfile(saveDir, [prefix ' - by vessel - AMLR strata.png']);
    print(ifile, '-dpng','-r300')
    crop_image(ifile)
   
    

    % Map coloured by stratum. Too crowded to be really useful...
    figure(3)
    clf
    plot_standard_map(strata, 'showStrataNames', false)

    % Use a different colour and symbol for each statum

    s = unique(nasc.Stratum);
    symbols = {'o' 'o' 'o' 'o' 'o' 'o' 'o' 'd' 'd' 'd' 'd' 'd' 'd' 'd' '<' '<'};
    h = nan(size(s));
    for i = 1:length(s)
        j = find(nasc.Stratum == s(i));
        h(i) = m_scatter(nasc.Longitude(j), nasc.Latitude(j), nasc.rho(j)/maxRho*maxSize+1, 'filled', symbols{i});
    end

    m_grid('box', 'on')
    legend(h, s, 'Location', 'SouthEast', 'NumColumns', 2, 'Interpreter', 'none')

    ifile = fullfile(saveDir, [prefix ' - by stratum.png']);
    print(ifile, '-dpng','-r300')
    crop_image(ifile)

    %%%%%%%%%%%
    % AMLR strata and transects
    figure(4)
    clf

    s = ["Bransfield" "Elephant" "Joinville" "West"];
    plot_standard_map(strata, 'centrePoint', [-58 -62], 'radius', 4, ...
        'strataToPlot', s, 'showStrataNames', true, ...
        'coastDetail', 'high', 'rotation', 0)

    for i = 1:length(s)
        j = find(nasc.Stratum == s(i));
        m_scatter(nasc.Longitude(j), nasc.Latitude(j), nasc.rho(j)/maxRho*maxSize+1, 'filled', 'o');

        % Transect label at northernmost point of each transect
%         t = unique(nasc.Transect(j));
%         for k = 1:length(t)
%             j = find(nasc.Stratum == s(i) & nasc.Transect == t(k));
%             [~, n_i] = max(nasc.Latitude(j));
%             m_text(nasc.Longitude(j(n_i)), nasc.Latitude(j(n_i)), nasc.Transect(j(n_i)));
%         end
    end
    plot_standard_map_rho_legend(legendScatterSizes, maxRho, maxSize)

    ifile = fullfile(saveDir, [prefix ' - AMLR.png']);
    print(ifile, '-dpng','-r300')
    crop_image(ifile)
    
    % and a version for the paper
    if strcmp(prefix, 'Krill density')
        ifile = fullfile(saveDir,  'Figure 8c.tiff');
        print(ifile, '-dtiff','-r1000')
        crop_image(ifile)
    end


    %%%%%%%%%%%
    % CCAMLR 2000 strata and transects
    figure(5)
    clf

    s = ["ESS" "Sand" "SG" "SS" "AP" "SSI" "SOI"];
    plot_standard_map(strata, 'centrePoint', [-45 -60], 'radius', 17.5, ...
        'strataToPlot', s, 'showStrataNames', true, ...
        'coastDetail', 'intermediate', 'showSubAreas', false, ...
        'rotation', 0)

    for i = 1:length(s)
        j = find(nasc.Stratum == s(i));
        m_scatter(nasc.Longitude(j), nasc.Latitude(j), nasc.rho(j)/maxRho*maxSize+1, 'filled', 'o');
    end
    plot_standard_map_rho_legend(legendScatterSizes, maxRho, maxSize)

    ifile = fullfile(saveDir,   [prefix ' - CCAMLR 2000.png']);
    print(ifile, '-dpng','-r300')
    crop_image(ifile)

    % and a version for the paper
    if strcmp(prefix, 'Krill density')
        ifile = fullfile(saveDir,  'Figure 8a.tiff');
        print(ifile, '-dtiff','-r1000')
        crop_image(ifile)
    end
    
    %%
    %%%%%%%%%%%
    % CCAMLR 2000 strata and transects as a density surface for comparison
    % to Figure 3 in Hewitt et. al., 2004 (DSRII)
    figure(6)
    clf

    s = ["ESS" "Sand" "SG" "SS" "AP" "SSI" "SOI"];
    %m_proj('Azimuthal Equal-area', 'long', -40, 'lat', -60, 'radius', 21, 'rectbox', 'on');
    m_proj('Transverse Mercator', 'long', [-64, -25], 'lat', [-51, -67], ...
        'clongitude', -40, 'rectboxt', 'on');
    patchColour = [0.5 0.5 0.5];
    m_gshhs_i('patch', 'FaceColor', patchColour);
    hold on

    % do a filled contour plot at these levels
    levels = [0, 5, 10, 20, 40, 80, 160, 1e10];
    % same colours are used by Hewitt
    colours = [254, 255, 167; % 0-5
        251, 235, 163; % 5-10
        247, 209, 151; % 10-20
        238, 166, 132; % 20-40
        226, 71, 96;  % 40-80
        204, 62, 95; % 80-160
        139, 67, 95]; % 160-infinity
    colours = colours / 255;
    colormap(colours)
    %colormap(jet(16))
    
    % Grid data for a filled contour plot
    j = contains(nasc.Stratum, s);
    
    grid_size = 0.1; % deg
    [X, Y] = meshgrid(min(nasc.Longitude(j)):grid_size:max(nasc.Longitude(j)), ...
                      max(nasc.Latitude(j)):-grid_size:min(nasc.Latitude(j)));
    F = scatteredInterpolant(nasc.Longitude(j), nasc.Latitude(j), nasc.rho(j), ...
        'natural', 'none');
    
    rho = F(X, Y);
    % then cut out regions outside the survey region (scatteredInterpolant
    % uses the convex hull which includes some areas we don't want to see).
    mask_strata = {'ESS', 'SS', 'AP'};
    mask = false(size(X)); % true means inside the strata
    for i = 1:length(strata.features)
        if sum(strcmp(strata.features(i).properties.stratum, mask_strata)) > 0
            poly = squeeze(strata.features(i).geometry.coordinates);
            in = inpolygon(X, Y, poly(:,1), poly(:,2));
            mask = mask | in;
        end
    end
    
    rho(~mask) = NaN;
    
    %m_contourf(X, Y, rho, levels, 'LineColor', 'none')
    % contourf somehow mucks up the colour bands, so do it more explicitly
    % using pcolor.
    C = NaN(size(X));
    for i = 1:length(levels)-1
        C(rho >= levels(i) & rho < levels(i+1)) = i;
    end
    
    h = m_image(X(1,:), Y(:,1), C);
    %set(h, 'EdgeColor', 'none');
    hold on
    m_plot(nasc.Longitude(j), nasc.Latitude(j), 'k.', 'MarkerSize', 5)
    m_grid('xtick', -70:10:-20, 'ytick', (-60:10:-50), 'box', 'on')
    %[ax, h] = m_contfbar([0.45, 0.93], 0.25, F(X, Y), levels, 'endpiece','yes', 'axfrac', 0.05);
    % and do the colourbar ticks manually cause it's the only way to get
    % what we want...
    h = colorbar('south');
    cb_lim = get(h, 'Limits');
    cb_step = (cb_lim(2)-1)/(length(levels)-1);
    set(h, 'Ticks', cb_lim(1) + (0:cb_lim(2)-1) * cb_step, 'TickLabels', levels(1:end-1))
    set(h, 'Position', [0.45 0.2 0.4 0.04])
    title(h, 'Krill density (g m^{-2})')
    
    ifile = fullfile(saveDir,   [prefix ' - CCAMLR 2000 - surface density.png']);
    print(ifile, '-dpng','-r300')
    crop_image(ifile)

    % and a version for the paper
    ifile = fullfile(saveDir,  'Figure 9b.tiff');
    print(ifile, '-dtiff','-r600')
    crop_image(ifile)

    %%
    %%%%%%%%%%%
    % Zoom on the South Orkney strata and transects
    figure(7)
    clf

    s = ["SOI" "SOC" "SOF"];
    plot_standard_map(strata, 'centrePoint', [-45.7 -60.75], 'radius', 2.5, ...
        'strataToPlot', s, 'showStrataNames', true, ...
        'coastDetail', 'fine', 'rotation', 0)

    for i = 1:length(s)
        j = find(nasc.Stratum == s(i));
        m_scatter(nasc.Longitude(j), nasc.Latitude(j), nasc.rho(j)/maxRho*maxSize+1, 'filled', 'o');
    end
    plot_standard_map_rho_legend(legendScatterSizes, maxRho, maxSize)

    ifile = fullfile(saveDir, [prefix ' - South Orkney.png']);
    print(ifile, '-dpng','-r300')
    crop_image(ifile)
    
    % and a version for the paper
    if strcmp(prefix, 'Krill density')
        ifile = fullfile(saveDir,  'Figure 8b.tiff');
        print(ifile, '-dtiff','-r1000')
        crop_image(ifile)
    end
    
    %%%%%%%%%%%%%
    % a map of where we are in the world
    figure(8)
    clf
    
    plot_standard_map(strata, 'showStrataNames', false, 'centrePoint', [-45 -60], ...
        'rotation', 0, 'radius', 30, 'strataColour', 'k')
    %m_coast('patch',[0.5 0.5 0.5]);
    m_grid('xticklabel',[],'yticklabel',[],'linestyle',':', 'box', 'off')
        
    ifile = fullfile(saveDir, 'globe with survey area.png');
    print(ifile, '-dpng','-r300')
    crop_image(ifile)

end


