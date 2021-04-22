
clear easypeasy
load(fullfile(resultsDir, 'SDWBA-TS-2019-038kHz.mat'), 'krill_ts')
easypeasy.length = [krill_ts.ts.ActualLength]'*1000;
easypeasy.TS038 = 10*log10([krill_ts.ts.sigma_avg]');

load(fullfile(resultsDir, 'SDWBA-TS-2019-070kHz.mat'), 'krill_ts')
easypeasy.TS070 = 10*log10([krill_ts.ts.sigma_avg]');

load(fullfile(resultsDir, 'SDWBA-TS-2019-120kHz.mat'), 'krill_ts')
easypeasy.TS120 = 10*log10([krill_ts.ts.sigma_avg]');

load(fullfile(resultsDir, 'SDWBA-TS-2019-200kHz.mat'), 'krill_ts')
easypeasy.TS200 = 10*log10([krill_ts.ts.sigma_avg]');


clear easypeasy05
load(fullfile(resultsDir, 'SDWBA-TS-2019-038kHz-0_5mm.mat'), 'krill_ts')
easypeasy05.length = [krill_ts.ts.ActualLength]'*1000;
easypeasy05.TS038 = 10*log10([krill_ts.ts.sigma_avg]');

load(fullfile(resultsDir, 'SDWBA-TS-2019-070kHz-0_5mm.mat'), 'krill_ts')
easypeasy05.TS070 = 10*log10([krill_ts.ts.sigma_avg]');

load(fullfile(resultsDir, 'SDWBA-TS-2019-120kHz-0_5mm.mat'), 'krill_ts')
easypeasy05.TS120 = 10*log10([krill_ts.ts.sigma_avg]');

load(fullfile(resultsDir, 'SDWBA-TS-2019-200kHz-0_5mm.mat'), 'krill_ts')
easypeasy05.TS200 = 10*log10([krill_ts.ts.sigma_avg]');


clf
plot(easypeasy.length, easypeasy.TS038, 'x-')
hold on
plot(easypeasy05.length, easypeasy05.TS038, '.-')

plot(easypeasy.length, easypeasy.TS070, 'x-')
plot(easypeasy05.length, easypeasy05.TS070, '.-')

plot(easypeasy.length, easypeasy.TS120, 'x-')
plot(easypeasy05.length, easypeasy05.TS120, '.-')

plot(easypeasy.length, easypeasy.TS200, 'x-')
plot(easypeasy05.length, easypeasy05.TS200, '.-')

legend('38','3805','70', '7005', '120','12005','200','20005')

TS038 = interp1(easypeasy.length, easypeasy.TS038, easypeasy05.length, 'pchip', 'extrap');
TS120 = interp1(easypeasy.length, easypeasy.TS120, easypeasy05.length, 'pchip', 'extrap');
TS200 = interp1(easypeasy.length, easypeasy.TS200, easypeasy05.length, 'pchip', 'extrap');

% Generate results in Stox TS table XML format...


