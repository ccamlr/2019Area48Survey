function day = dayOrNight(lon, lat, date)
% Tell if the given datetime & position is in civil day or civil night. 
%     lon and lat should be in decimal degrees (-ve for S and W)
%     date is always taken to be in UTC and can be in any format accepted
%     by datenum()
% 
% Adapted from sunRiseSet.m by Richard Droste
% https://github.com/rdroste/sunRiseSet
% 
% Copyright (c) 2017, Richard Droste
% All rights reserved.
% 
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are
% met:
% 
%     * Redistributions of source code must retain the above copyright
%       notice, this list of conditions and the following disclaimer.
%     * Redistributions in binary form must reproduce the above copyright
%       notice, this list of conditions and the following disclaimer in
%       the documentation and/or other materials provided with the distribution
% 
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
% POSSIBILITY OF SUCH DAMAGE.

% Reverse engineered from the NOAA Excel spreadsheet:
% (https://www.esrl.noaa.gov/gmd/grad/solcalc/calcdetails.html)
% 
% The formulas are from:
% Meeus, Jean H. Astronomical algorithms. Willmann-Bell, Incorporated, 1991.

nDays = floor(datenum(date)-datenum('30-dec-1899'));

nTimes = 24*3600;
tArray = seconds(timeofday(datetime(date)))/nTimes;

UTCoff = 0; % Fixed to UTC timezone

% Compute
% Letters correspond to colums in the NOAA Excel
E = tArray;
F = nDays+2415018.5+E-UTCoff/24;
G = (F-2451545)/36525;
I = mod(280.46646+G.*(36000.76983+G*0.0003032),360);
J = 357.52911+G.*(35999.05029-0.0001537*G);
K = 0.016708634-G.*(0.000042037+0.0000001267*G);
L = sind(J).*(1.914602-G.*(0.004817+0.000014*G))+sind(2*J).* ...
    (0.019993-0.000101*G)+sind(3*J)*0.000289;
M = I+L;
P = M-0.00569-0.00478*sind(125.04-1934.136*G);
Q = 23+(26+((21.448-G.*(46.815+G.*(0.00059-G*0.001813))))/60)/60;
R = Q+0.00256*cosd(125.04-1934.136*G);
T = asind(sind(R).*sind(P));
U = tand(R/2).*tand(R/2);
V = 4*rad2deg(U.*sin(2*deg2rad(I))-2*K.*sin(deg2rad(J))+4*K.*U.*sin(deg2rad(J)).* ...
    cos(2*deg2rad(I))-0.5.*U.*U.*sin(4*deg2rad(I))-1.25.*K.*K.*sin(2.*deg2rad(J)));
AB = mod(E*1440+V+4*lon-60*UTCoff,1440);
if AB/4 < 0, AC = AB/4+180;else, AC = AB/4-180; end
AD = acosd(sind(lat)*sind(T)+cosd(lat)*cosd(T).*cosd(AC));
W = acosd(cosd(90.833)./(cosd(lat)*cosd(T))-tand(lat)*tand(T));
X = (720-4*lon-V+UTCoff*60)*60;

elev_ang_corr = 90-AD;

if elev_ang_corr > -6 % the defintion of civil time.
    day = true;
else
    day = false;
end




    
