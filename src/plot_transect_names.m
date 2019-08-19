function plot_transect_names(nasc, s)
    
for i = 1:length(s)
    j = find(nasc.Stratum == s(i));
    m_plot(nasc.Longitude(j), nasc.Latitude(j), 'k.');
    hold on
    
    % and a transect label at the northernmost part
    t = unique(nasc.Transect(j));
    for k = 1:length(t)
        j = find(nasc.Stratum == s(i) & nasc.Transect == t(k));
        [~, n_i] = max(nasc.Latitude(j));
        m_text(nasc.Longitude(j(n_i)), nasc.Latitude(j(n_i)), nasc.Transect(j(n_i)), ...
            'FontSize', 8, 'Color', [ 0.2422    0.1504    0.6603], ...
            'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center');
    end
end

