%%

% make some plots that compare the krill density from integrations done
% with a -30 and a -40 dB threshold in the CCMALR echoview template

do_define_directories

load(fullfile(resultsDir2020, 'Final results - swarm'), 'results');
results2019 = load(fullfile(resultsDir2019, 'Final results - swarm'), 'results');
results2019 = results2019.results;

[N2019, edges2019] = histcounts(log10(results2019.nasc.rho));
[N2020, edges2020] = histcounts(log10(results.nasc.rho));
distance2019 = (1:length(results2019.nasc.rho)) * 1.852;
distance2020 = (1:length(results.nasc.rho)) * 1.852;

clf
% subplot(2,1,1)
% plot(edges2019(1:end-1), N2019, 'k', 'LineWidth', 1.5)
% hold on
% plot(edges2020(1:end-1), N2020, 'k:', 'LineWidth', 1.5)
% set(gca, 'XTickLabel', [0.001 0.01 0.1 1 10 100 1000 10000])
% %title('Krill density distribution')
% xlabel('Krill density (g m^{-2})')
% textLoc('(a)', 'NorthWest');
% ylabel('Proportion')
% %legend('-40 dB', '-30 dB', 'Location', 'NorthWest')

%subplot(2,1,2)
plot(distance2019, cumsum(results2019.nasc.rho), 'k:', 'LineWidth', 1.5)
hold on
plot(distance2020, cumsum(results.nasc.rho), 'k', 'LineWidth', 1.5)
%title('Cumulative krill density')
%legend('-40 dB', '-30 dB', 'Location', 'NorthWest')
ylabel(['Cumulative krill density (g m^{' char(8211) '2})'])
xlabel('Cumulative transect distance (km)')
%textLoc('(b)', 'NorthWest');

ifile = fullfile(resultsDir, 'Threshold rho comparison.png');
print(ifile, '-dpng','-r300')
crop_image(ifile)
    
% and a version for the paper
ifile = fullfile(resultsDir, 'Figure 7.tiff');
print(ifile, '-dtiff','-r1000')
crop_image(ifile)

    