function do_location_maps(strata, saveDir)
% Make various plots of the survey area. 


    %%%%%%%%%%%
    % AMLR strata and transects
    figure(1)
    clf

    s = ["Bransfield" "Elephant" "Joinville" "West"];
    plot_standard_map(strata, 'centrePoint', [-58 -62], 'radius', 4, ...
        'strataToPlot', s, 'showStrataNames', true, ...
        'coastDetail', 'high', 'rotation', 0)
    m_grid('box', 'on')
    
    ifile = fullfile(saveDir,  'Figure 1d.png');
    print(ifile, '-dpng','-r600')
    crop_image(ifile)

    %%%%%%%%%%%
    % CCAMLR 2000 strata and transects
    figure(2)
    clf

    s = ["ESS" "Sand" "SG" "SS" "AP" "SSI" "SOI"];
    plot_standard_map(strata, 'centrePoint', [-45 -60], 'radius', 17.5, ...
        'strataToPlot', s, 'showStrataNames', true, ...
        'coastDetail', 'intermediate', 'showSubAreas', true, ...
        'rotation', 0)
    m_grid('box', 'on')
    
    ifile = fullfile(saveDir,  'Figure 1b.png');
    print(ifile, '-dpng','-r600')
    crop_image(ifile)
    
    %%%%%%%%%%%
    % Zoom on the South Orkney strata and transects
    figure(3)
    clf

    s = ["SOI" "SOC" "SOF"];
    plot_standard_map(strata, 'centrePoint', [-45.7 -60.75], 'radius', 2.5, ...
        'strataToPlot', s, 'showStrataNames', true, ...
        'coastDetail', 'fine', 'rotation', 0)
    m_grid('box', 'on')

    ifile = fullfile(saveDir,  'Figure 1c.png');
    print(ifile, '-dpng','-r600')
    crop_image(ifile)
    
    %%%%%%%%%%%%%
    % a map of where we are in the world
    figure(4)
    clf
    
    plot_standard_map(strata, 'showStrataNames', false, 'centrePoint', [-45 -60], ...
        'rotation', 0, 'radius', 30, 'strataColour', 'k')
    %m_coast('patch',[0.5 0.5 0.5]);
    m_grid('xticklabel',[],'yticklabel',[],'linestyle',':', 'box', 'off')
        
    ifile = fullfile(saveDir,  'Figure 1a.png');
    print(ifile, '-dpng','-r600')
    crop_image(ifile)
    
end


