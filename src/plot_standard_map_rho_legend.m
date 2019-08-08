function plot_standard_map_rho_legend(legendScatterSizes, maxRho, maxSize)


legentry = cell(size(legendScatterSizes));
bubleg = nan(size(legendScatterSizes));
for ind = 1:numel(legendScatterSizes)
   bubleg(ind) = m_plot(-45, -60, 'ko', 'MarkerSize', sqrt(legendScatterSizes(ind)/maxRho*maxSize+1), ...
       'MarkerFaceColor', 'k');
   hold on
   set(bubleg(ind),'visible','off')
   legentry{ind} = num2str(legendScatterSizes(ind));
end

m_grid('box', 'on')

leg = legend(bubleg, legentry, 'Location', 'SouthEast');
title(leg, '\rho (g\cdotm^{-2})')
set(leg, 'fontWeight', 'normal')
legend('boxoff')
