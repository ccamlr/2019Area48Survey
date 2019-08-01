% calculate the C for conversion of NASC to areal density
clear all; clf

% load krill length-frequencies by cluster
load clusters123_sizesCCAMLR_1.txt 
lf=clusters123_sizesCCAMLR_1';

c=1456;                             % sound speed
freq=[38 120 200]*1e3;              % acoustic frequency
k=(2*pi.*freq)/c;                   % acoustic wavenumber
lo=38.35;                           % reference length (mm)
Lo=38.35*1e-3;                      % reference length (m)
l=lf(1,:);                          % krill lengths (mm)
L=l*1e-3;                           % krill lengths (m)
w=2.236e-3*l.^3.314;                % krill mass (g) (CCAMLR 2000)

% load SDWBA TS strengths for krill lengths
load TS_krill_length_values_alt_fin.mat

% compute and plot SDWBA TS model
for i=1:size(freq,2),
    kL(i,:)=k(i)*L;
    LN(i,:)=L;
    wN(i,:)=w;
end

wNw=wN.*lf(2:4,:);  % frequency-weighted weights

% Set TS from loaded file of TS krill lengths
Fullm_TS = T_TS(1:3,:);

hold on
col=['m' 'b' 'r' 'g'];

fig1=figure(1);
set(fig1,'Color',[1 1 1]);

for cluster=1:1:3,
    for i=1:size(freq,2),
        SDWBA_full(cluster,i)=sum(lf(cluster+1,:).*wN(i,:))/sum(lf(cluster+1,:).*(4*pi*10.^(Fullm_TS(i,:)/10)));
    end
end
SDWBA_full=SDWBA_full'/1e+3/1852^2          % Convert C to units of m^2/n.mi.^2

% Calculate Greene et al. model for comparison
GTS(1,:)=-132.44+34.85*log10(l);
GTS(2,:)=-127.45+34.85*log10(l);
GTS(3,:)=-125.23+34.85*log10(l);

for cluster=1:1:3,
    for i=1:size(freq,2),
        Greene_CF(cluster,i)=sum(lf(cluster+1,:).*wN(i,:))/sum(lf(cluster+1,:).*(4*pi*10.^(GTS(i,:)/10)));
    end
end
Greene_CF=Greene_CF'/1e+3/1852^2	% Convert C to units of m^2/n.mi.^2

line(l,Fullm_TS(2,:),'LineWidth',2,'LineStyle','-','Color',col(3))
line(l,GTS(2,:),'LineWidth',2,'LineStyle','-','Color',col(4))

xlabel('Krill Length (mm)')
ylabel('Target Strength (dB)')
grid

var=get(fig1,'children');
set(var,'xscale','log');
axis([10,100,-95,-65])
ylabel('TS (dB)')
xlabel('Length (mm)')
legend('SDWBA full','Greene \it{et al.}',0)
grid on
box on
