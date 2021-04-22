%%

figure(1)
clf
strata.features = [];

m_proj('Azimuthal Equal-area', ...
    'long', -40, ...
    'lat', -59, ...
    'radius', 12, ...
    'rectbox', 'on');

m_gshhs_i
hold on

d = readtable('../map_data/APOP stations.csv');
i = str2double(d.toktnr) == 2019701;
j = str2double(d.toktnr) == 2019828;

c1 = [216,179,101]/255;
c2 = [90,180,172]/255;
m_plot(str2double(d.lengde(i)), str2double(d.bredde(i)), 'o', 'Color', c1, 'MarkerFaceColor', c1)

m_plot(str2double(d.lengde(j)), str2double(d.bredde(j)), '^', 'Color', c2, 'MarkerFaceColor', c2)

m_grid('box', 'on')
 
ifile='APOP_positions.png';
print(ifile, '-dpng','-r300')
crop_image(ifile)
