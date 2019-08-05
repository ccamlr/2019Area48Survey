function plot_standard_map(strata, varargin)
% Draw an appropriate map background for the survey, using the m_map
% toolbox. 

p = inputParser;
p.addOptional('showStrataNames', true)
parse(p, varargin{:});

m_proj('Azimuthal Equal-area', 'lat', -60, 'long', -45, 'radius', 17, ...
    'rectbox', 'on');
m_gshhs_l('patch', [.5 .5 .5]);
hold on

for i = 1:length(strata.features)
    poly = squeeze(strata.features(i).geometry.coordinates);
    m_line(poly(:,1), poly(:,2), 'color', [0.5 0.5 0.5])
    
    if p.Results.showStrataNames
        % Add the stratum name at the most northly point of the stratum
        [~, j] = max(poly(:,2));
        m_text(poly(j(1),1), poly(j(1),2), strata.features(i).properties.stratum, ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom')
    end
end


