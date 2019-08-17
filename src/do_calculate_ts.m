%%
% Calculate the krill TS using parameters relevant to the 2019 survey. 
%
% Only need to run this code once - the results are saved and then used as
% necessary.

% Storage parameters
do_define_directories
SDWBAdir = fullfile(baseDir, repoDir, 'src', 'SDWBApackage2010');

resultsFile = 'SDWBA-TS-2019-120kHz';

addpath(SDWBAdir)

% Key parameters
ActualLengths = (10:67)*1e-3; % [m]
h0 = 1.0279; % sound speed contrast
g0 = 1.0357; % density contrast
c = 1456; % sound speed [m/s]

% Operational parameters
frequency = 120e3; % [Hz]
phi = -90:1:269; % [deg]
noise_realisations = 100;

% Basic parameters
Animal_shape_file = 'Esuperba_AT38_35mm_McGeheeetal_1998';
load(fullfile(SDWBAdir, Animal_shape_file), 'a', 'r')
a0 = a;
r0 = r;
N0 = length(a)-1;
L0 = 38.35*1e-3; % [m]
fatness_factor = 1.4;
freq0 = 120e3; % [Hz]
stdphase0 = sqrt(2)/2;

% For tilt averaging
meanorientation = -20; % [deg]
stdorientation = 28; % [deg]

% Do the TS calculations...
phi = 90 + phi; 
sigma_avg = nan(size(ActualLengths));

% Store the overall parameters
clear results
results.parameters = struct('Animal_shape_file', Animal_shape_file, 'Processing_Date', datestr(now), ...
    'L0', L0, 'fatness_factor', fatness_factor, 'N0', N0, 'stdphase0', stdphase0, ...
    'freq0', freq0, 'g0', g0, 'h0', h0, ...
    'noise_realisations', noise_realisations, 'frequency', frequency, ...
    'phi', phi, 'c', c);

% Setup a plot that is updated with each new TS point
figure(1)
clf
xlabel('Length (mm)')
ylabel('TS (dB re 1m^2)')
title(['Krill backscatter at ' num2str(results.parameters.frequency*1e-3) ' kHz'])
h = plot(NaN, NaN, '.-');

% Calculate the TS at each length
for i_len = 1:length(ActualLengths)
    disp(['Calculating length ' num2str(ActualLengths(i_len)*1e3) ' mm.'])

    kL = 2*pi*frequency / c * ActualLengths(i_len) ;        % wave number
    
    scaling_factor = ActualLengths(i_len)/L0 ;
    r = r0 * scaling_factor ;
    a = a0 * scaling_factor * fatness_factor;
   
    % Allocate a variable to accumulate sigma_bs into
    BSsigmatot = zeros([size(phi),noise_realisations]);

    % Do multiple model runs
    for irealisation = 1:noise_realisations
        %disp(['  Calculating ' num2str(irealisation) ' of ' num2str(noise_realisations)])
        [~, BSsigma, ~, Stdphase_vs_freq, Cylinders_vs_freq] = BSTS_SDWBA_2010(frequency,r,a,h0,g0,phi,[stdphase0 freq0 N0],c,[],'');
        BSsigmatot(:,:,irealisation) = BSsigma;
    end
    
    BSsigma_StandardDeviation = std(BSsigmatot, [], 3);
    BSsigma = mean(BSsigmatot, 3);

    % work out average TS
    orientation = GaussianOrientation(phi, 90-meanorientation, stdorientation);
    [sigma, ~] = AverageTSorientation(BSsigma, orientation, phi);
    
    % Store results in a structure
    results.ts(i_len) = struct('scaling_factor', scaling_factor, ...
        'BSsigma', BSsigma, ...
        'BSsigma_StandardDeviation', BSsigma_StandardDeviation, ...
        'ActualLength', ActualLengths(i_len), ...
        'a', a, 'r', r, ...
        'Stdphase_vs_freq', Stdphase_vs_freq, ...
        'Cylinders_vs_freq', Cylinders_vs_freq, ...
        'sigma_avg', sigma);
    
    % save each iteration in case of crash
    krill_ts = results;
    save(fullfile(resultsDir, resultsFile), 'krill_ts')

    % update the plot. Doing it this way doesn't give focus to the plot
    set(h, 'XData', [results.ts.ActualLength]*1e3, ...
        'YData', 10*log10([results.ts.sigma_avg]))
    drawnow()

end

rmpath(SDWBAdir)

%%
% and a little comparison between what is calculated above and what is in
% Table 2 of WG-EMM-16/38 (and is also in the file
% ASAM_TS_krill_length_values_alt_fin.mat, which is stored in the results
% directory)

figure(2)
clf

do_define_directories

clear ts
%ts = readtable(fullfile(resultsDir, 'sdwba_ts.csv'));
load(fullfile(resultsDir, 'ASAM_TS_krill_length_values_alt_fin'), 'T_TS', 'krill_ls')
% and make the data structure a little easier to understand...
ts.ts038 = T_TS(1,:)';
ts.ts120 = T_TS(2,:)';
ts.ts200 = T_TS(3,:)';
ts.length = krill_ls';
% That dataset goes to 65 mm, but we have krill to 67mm, so deal with that
% in the plots
l = 10:67;

load(fullfile(resultsDir, 'SDWBA-TS-2019-038kHz.mat'), 'krill_ts')
plot(l, 10*log10([krill_ts.ts.sigma_avg]') - [ts.ts038; NaN; NaN])
hold on

load(fullfile(resultsDir, 'SDWBA-TS-2019-120kHz.mat'), 'krill_ts')
plot(l, 10*log10([krill_ts.ts.sigma_avg]') - [ts.ts120; NaN; NaN])

load(fullfile(resultsDir, 'SDWBA-TS-2019-200kHz.mat'), 'krill_ts')
plot(l, 10*log10([krill_ts.ts.sigma_avg]') - [ts.ts200; NaN; NaN])

xlabel('Length (mm)')
ylabel('\Delta TS of krill (dB re 1m^2)')
title({'Difference between WG-EMM-16/38 (Table 2)', ' and SDWBApackage2010'})
grid
legend('38 kHz','120 kHz','200 kHz')
legend boxoff

print(fullfile(resultsDir, 'SDWBA_TS_comparison'), '-dpng', '-r300')
