function [rho_j, L_j] = calcDensity(data)
   % Calculates mean areal krill biomass density for a transect as per
   % equation 9 of EMM-16/38
   %
   % Input:
   %  Table with variables:
   %    - Latitude [decimal degrees]
   %    - Longitude [decimal degrees]
   %    - NASC (Nautical area scattering coefficient, also s_A [m^2 / nmi^2])
   %    - C (conversion factor [g/m^2])
   %
   % All rows in the supplied table are used, and assumed to be from the
   % one transect.
   %
   % The NASC values are weighted to account for distance off the nominal
   % straight-line transect with the assumption that the transects are more
   % or less N/S transects, rather than E-W transects. 
   %
   % Output:
   %   - rho_j (mean areal krill biomass density [g/m^2])
   %   - L_j (length of the transect [nmi])
   

   % Find extent of transect in N/S direction
   [~, i1] = max(data.Latitude);
   [~, i2] = min(data.Latitude);
   
   t = struct('lat1', data.Latitude(i1), 'lon1', data.Longitude(i1), ...
       'lat2', data.Latitude(i2), 'lon2', data.Longitude(i2));

   % Find the difference between the start and end latitude of the transect
   delta_initial = t.lat1 - t.lat2;
    
   % Find the distance in nautical miles between start and end points using
   % the Haversine formula
   [~, nmi, ~, transect_course] = haversine([t.lat1 t.lon1], [t.lat2 t.lon2]);
    
   % Compute expected change in latitude for each nautical mile
   e_d_lat = delta_initial/nmi ;
    
   % Total # of intervals
   numInts = height(data);
    
   W_I = nan(numInts, 1); % Interval lengths (n.mi.)
    
   % Compute change in latitude from mile to mile. We can't guarantee that
   % the data are in spatial order (since some transects were done by more
   % than one vessel), so sort on latitude first
   [s_lat, i] = sort(data.Latitude);
   s_lon = data.Longitude(i);
   
   d_lat = abs(diff(s_lat));
   
   % Calculate weighting factor (Equation 10)
   W_I = (abs(e_d_lat) - abs(e_d_lat-d_lat)) / abs(e_d_lat); % [1]
  
   % If deviation from standard track line is < 10%, then simply set
   % weighting factor to 1
   W_I(W_I >= 0.9) = 1;
   
   % and an addition that wasn't in the CCAMLR 2000 procedure for when W_I
   % is negative. The 2000 procedure used start and end latitudes for each
   % 1 n.mile long NASC bin that were provided by Echoview. However, 
   % the exports we get from the current Echoview
   % template don't provide this - we only get one latitude per nasc bin
   % (probably the middle point.
   % For the current data, the actual latitude change in a nasc bin is
   % calculated from 
   % the difference in latitude to the next nasc bin. But, when there are
   % gaps in the transects (due to not exactly resuming the transect after
   % a station, for example), the difference in latitude can be much greater
   % than the expected and this leads to W_I's that are negative. The use of
   % W_I for weighting assumes that the W_I's are always between 0 and 1,
   % so this is a problem.
   % The solution is:
   % - assume that transects are straight lines in the lat/log coordinate
   %   system
   % - calculate the expected direction for the transect
   % - for each nasc value estimate the direction to the next nasc value
   % - use the difference in directions to estimate W_I's. 
   % - this should give near identical results to the CCAMLR 2000 procedure
   %   if using start and stop positons for each nasc value.
   for i = 1:height(data)-1
       [~, ~, ~, c] = haversine([s_lat(i) s_lon(i)], ...
           [s_lat(i+1) s_lon(i+1)]);
       course(i) = c;
   end
   
   W_I = ones(numInts, 1)
   
   % Compute transect length (n.mi.) by summing interval lengths
   L_j = sum(W_I);
    
   % Convert NASC to krill biomass density using weighting and conversion
   % factors. 
   % This is equation 9 in EMM-16/38, corrected by using s_a instead of s_A
   % (which involves dividing by 1852^2). The 4pi is accounted for by using
   % sigma instead of sigma_bs in C.
   rho_j = sum(data.NASC / 1852^2.* data.C .* W_I) ./ L_j;
end

