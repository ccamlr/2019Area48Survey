function plot_standard_map(strata, varargin)
% Draw an appropriate map background for the survey, using the m_map
% toolbox. 

p = inputParser;
p.addOptional('showStrataNames', true)
p.addOptional('centrePoint', [-45 -60])
p.addOptional('radius', 17)
p.addOptional('strataToPlot', []);
p.addOptional('coastDetail', 'low');
p.addOptional('showSubAreas', false)
p.addOptional('strataColour', [0.5 0.5 0.5])

parse(p, varargin{:});

m_proj('Azimuthal Equal-area', ...
    'long', p.Results.centrePoint(1), ...
    'lat', p.Results.centrePoint(2), ...
    'radius', p.Results.radius, ...
    'rectbox', 'on');

patchColour = [0.5 0.5 0.5];
switch p.Results.coastDetail
    case 'coarse'
        m_gshhs_c('patch', patchColour);
    case 'low'
        m_gshhs_l('patch', patchColour);
    case 'intermediate'
        m_gshhs_i('patch', patchColour);
    case 'high'
        m_gshhs_h('patch', patchColour);
    case 'fine'
        m_gshhs_f('patch', patchColour);
end

hold on

if p.Results.showSubAreas
    b = m_shaperead('../map_data/asd-shapefile-WGS84');
    % only want to show specific sub areas
    areas = {'Subarea 48.1' 'Subarea 48.2' 'Subarea 48.3' 'Subarea 48.4' };
    for i = 1:length(b.ncst)
        if ismember(b.dbfdata{i,2}, areas)
            m_line(b.ncst{i}(:,1), b.ncst{i}(:,2), 'color', 'k')
            m_text(min(b.ncst{i}(:,1)), max(b.ncst{i}(:,2)), b.dbfdata{i,2}(9:end), ...
                'HorizontalAlignment', 'left', 'VerticalAlignment', 'top')
        end
    end
end

for i = 1:length(strata.features)
    if isempty(p.Results.strataToPlot) || ...
        contains(strata.features(i).properties.stratum,p.Results.strataToPlot)
        poly = squeeze(strata.features(i).geometry.coordinates);
        m_line(poly(:,1), poly(:,2), 'color', p.Results.strataColour)
    
        if p.Results.showStrataNames
            % Add the stratum name at the most northly point of the stratum
            [~, j] = max(poly(:,2));
            m_text(poly(j(1),1), poly(j(1),2), strata.features(i).properties.stratum, ...
                'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom')
        end
    end
end


