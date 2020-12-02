%%

% make some plots that compare the krill density from integrations done
% with a -30 and a -40 dB threshold in the CCMALR echoview template

do_define_directories

load(fullfile(resultsDir, 'Final results - swarm'), 'results');
results2019 = load(fullfile(resultsDir2, 'Final results - swarm'), 'results');
results2019 = results2019.results;

[N2019, edges2019] = histcounts(log10(results2019.nasc.rho));
[N2020, edges2020] = histcounts(log10(results.nasc.rho));

clf
subplot(2,1,1)
plot(edges2019(1:end-1), N2019)
hold on
plot(edges2020(1:end-1), N2020)
title('Krill density distribution')
xlabel('log_{10}(krill density) (g m^{-2})')
ylabel('Counts')
legend('-40 dB', '-30 dB', 'Location', 'NorthWest')

subplot(2,1,2)
plot(cumsum(results2019.nasc.rho))
hold on
plot(cumsum(results.nasc.rho))
title('Cumulative krill density')
legend('-40 dB', '-30 dB', 'Location', 'NorthWest')
ylabel('Krill density (g m^{-2})')
xlabel('NASC number')

ifile = fullfile(resultsDir, 'Threshold rho comparison.png');
print(ifile, '-dpng','-r300')
crop_image(ifile)
    
    