function do_day_night_coverage_maps(nasc_all, nasc_day, strata, prefix, saveDir)
% Make various plots of the krill areal density. 

    sizeAll = 5; % [points^2] of drawn circles
    sizeDay = 10; % [points^2] of drawn circles

    %%%%%%%%%%%
    % AMLR strata and transects
    figure(1)
    clf

    s = ["Bransfield" "Elephant" "Joinville" "West"];
    plot_standard_map(strata, 'centrePoint', [-58 -62], 'radius', 4, ...
        'strataToPlot', s, 'showStrataNames', false, ...
        'coastDetail', 'high', 'rotation', 0)

    for i = 1:length(s)
        j = find(nasc_all.Stratum == s(i));
        m_scatter(nasc_all.Longitude(j), nasc_all.Latitude(j), sizeAll, 'filled', 'ko');
        
        j = find(nasc_day.Stratum == s(i));
        m_scatter(nasc_day.Longitude(j), nasc_day.Latitude(j), sizeDay, 'filled', 'go');    
    end
    m_grid('box', 'on')
    
    ifile = fullfile(saveDir,  'Figure 5b.tiff');
    print(ifile, '-dtiff','-r1000')
    crop_image(ifile)

    %%%%%%%%%%%
    % CCAMLR 2000 strata and transects
    figure(2)
    clf

    s = ["ESS" "Sand" "SG" "SS" "AP" "SSI" "SOI"];
    plot_standard_map(strata, 'centrePoint', [-45 -60], 'radius', 17.5, ...
        'strataToPlot', s, 'showStrataNames', false, ...
        'coastDetail', 'intermediate', 'showSubAreas', false, ...
        'rotation', 0)

    for i = 1:length(s)
        j = find(nasc_all.Stratum == s(i));
        m_scatter(nasc_all.Longitude(j), nasc_all.Latitude(j), sizeAll, 'filled', 'ko');

        j = find(nasc_day.Stratum == s(i));
        m_scatter(nasc_day.Longitude(j), nasc_day.Latitude(j), sizeDay, 'filled', 'go');
    end
    m_grid('box', 'on')
     
    ifile = fullfile(saveDir,   'Figure 5a.tiff');
    print(ifile, '-dtiff','-r1000')
    crop_image(ifile)

    %%%%%%%%%%%
    % Zoom on the South Orkney strata and transects
    figure(3)
    clf

    s = ["SOI" "SOC" "SOF"];
    plot_standard_map(strata, 'centrePoint', [-45.7 -60.75], 'radius', 2.5, ...
        'strataToPlot', s, 'showStrataNames', false, ...
        'coastDetail', 'fine', 'rotation', 0)

    for i = 1:length(s)
        j = find(nasc_all.Stratum == s(i));
        m_scatter(nasc_all.Longitude(j), nasc_all.Latitude(j), sizeAll, 'filled', 'ko');

        j = find(nasc_day.Stratum == s(i));
        m_scatter(nasc_day.Longitude(j), nasc_day.Latitude(j), sizeDay, 'filled', 'go');
    end
    m_grid('box', 'on')
     
    ifile = fullfile(saveDir, 'Figure 5c.tiff');
    print(ifile, '-dtiff','-r1000')
    crop_image(ifile)
 
end


