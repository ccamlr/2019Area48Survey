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
    
    %%%%%%%%%%%
    % Zoom on the South Orkney strata and transects
    figure(6)
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
    figure(7)
    clf
    
    plot_standard_map(strata, 'showStrataNames', false, 'centrePoint', [-45 -60], ...
        'rotation', 0, 'radius', 30, 'strataColour', 'k')
    %m_coast('patch',[0.5 0.5 0.5]);
    m_grid('xticklabel',[],'yticklabel',[],'linestyle',':', 'box', 'off')
        
    ifile = fullfile(saveDir, 'globe with survey area.png');
    print(ifile, '-dpng','-r300')
    crop_image(ifile)

end


