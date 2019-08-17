function [krill_len, krill_sigma] = get_krill_sigma(resultsDir)
% Return krill sigma and length at 120 kHz

load(fullfile(resultsDir, 'SDWBA-TS-2019'), 'krill_ts')
krill_len = [krill_ts.ts(:).ActualLength]; % [m]
krill_sigma = [krill_ts.ts(:).sigma_avg]*4*pi; % [m^2] Why oh why are the CCAMLR equations not using backscattered sigma????

