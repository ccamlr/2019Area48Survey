%%

% Pick out the water properties to use for the calibration, from the
% seabird CTD cast

baseDir = 'I:\KRILL2019\cruise_folders\KJH2019999';
d = ctd_rd(fullfile(baseDir, 'PHYSICS\CTD\CTD data/CAL_SBE37SMP-RS232_03720469_2019_03_07.cnv'),'LAB');
clf
s = d.sal00;
t = d.tv290C;
d = d.depSM;

[~, j] = max(d); % max depth
i = find(d > 1); % deeper than 1m
i = i(i >= j); % on the upcast

d = d(i);
s = s(i);
t = t(i);

subplot(1,2,1)
plot(t, d,'.-')
set(gca,'Ydir','reverse')
xlabel('Temperature (m/s)')
ylabel('Depth (m)')

subplot(1,2,2)
plot(s, d,'.-')
set(gca,'Ydir','reverse')
xlabel('Salinity (PSU)')
ylabel('Depth (m)')

% cal sphere was at about 21 m deep. Transducer was 5 m deep.
calDepth = 21; % [m]

[~, i] = min(abs(d - calDepth));

t(i)
s(i)

c = sw_svel(s(i), t(i), calDepth);
alpha38 = sw_absorption(38, s(i), t(i), calDepth, 'fandg', 8.0)*1e-3;
alpha120 = sw_absorption(120, s(i), t(i), calDepth, 'fandg', 8.0)*1e-3;

% c = 1452.6;
% t(i) = 0.9348
% s(i) = 34.1578

% The calibration was done with the environment variables to set 10degC, 24
% PSU, so need to change these to the values obtained above

rawDir = fullfile(baseDir, 'ACOUSTIC\EK60\EK60_CALIBRATION\20190307');
modDir = fullfile(baseDir, 'ACOUSTIC\EK60\EK60_CALIBRATION\20190307_modified');

d = dir(fullfile(rawDir, '*.raw'));
for i = 1:length(d)
    [header, data] = readEKRaw(fullfile(d(i).folder, d(i).name));
    data.pings(1).soundvelocity = c * ones(size(data.pings(1).soundvelocity));
    data.pings(1).absorptioncoefficient = alpha38 * ones(size(data.pings(1).absorptioncoefficient));
    data.pings(2).soundvelocity = c * ones(size(data.pings(2).soundvelocity));
    data.pings(2).absorptioncoefficient = alpha120 * ones(size(data.pings(2).absorptioncoefficient));
    writeEKRaw(fullfile(modDir, d(i).name), header, data);    
end

