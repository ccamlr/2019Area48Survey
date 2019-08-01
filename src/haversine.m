function [km, nmi, mi, course] = haversine(loc1, loc2)
    % HAVERSINE Compute distance between locations using Haversine formula
    % KM = HAVERSINE(LOC1, LOC2) returns the distance KM in km between
    % locations LOC1 and LOC2 using the Haversine formula. LOC1 and LOC2 are
    % latitude and longitude coordinates that can be expressed as decimal
    % degrees (where negative indicates West/South).
    %
    % [KM, NMI, MI] = HAVERSINE(LOC1, LOC2) returns the computed distance in
    % kilometers (KM), nautical miles (NMI), and miles (MI).
    %
    %
    % The first element indicates the latitude while the second is the
    % longitude.
    %
    % Notes
    % The Haversine formula is used to calculate the great-circle
    % distance between two points, which is the shortest distance over
    % the earth's surface.
    %
    % This program was created using equations found on the website
    % http://www.movable-type.co.uk/scripts/latlong.html
    % Created by Josiah Renfree
    % May 27, 2010
    % 
    % Modified by Gavin Macaulay, August 2019 to provide course and to
    % remove the parsing of lat/lon (only accepts decimal degrees now).
    
    % Check user inputs
    % If two inputs are given, display error
    if ~isequal(nargin, 2)
        error('User must supply two location inputs')
        % If two inputs are given, handle data
    else
        locs = {loc1 loc2}; % Combine inputs to make checking easier
    end
    
    % Check that both cells are a 2-valued array
    if any(cellfun(@(x) ~isequal(length(x),2), locs))
        error('Incorrect number of input coordinates')
    end
    
    % Convert all decimal degrees to radians
    locs = cellfun(@(x) x .* pi./180, locs, 'UniformOutput', 0);
    
    % Begin calculation
    R = 6371; % Earth's radius in km
    delta_lat = locs{2}(1) - locs{1}(1); % difference in latitude
    delta_lon = locs{2}(2) - locs{1}(2); % difference in longitude
    a = sin(delta_lat/2)^2 + cos(locs{1}(1)) * cos(locs{2}(1)) * ...
        sin(delta_lon/2)^2;
    c = 2 * atan2(sqrt(a), sqrt(1-a));
    km = R * c; % distance in km
    % Convert result to nautical miles and miles
    nmi = km * 0.539956803; % nautical miles
    mi = km * 0.621371192; % miles
    
    % and the bearing bearing at the start position that will get you to
    % the end position if travelling along a great circle line
    theta = atan2(sin(delta_lon)*cos(locs{2}(1)), ...
        cos(locs{1}(1))*sin(locs{2}(1))-...
        sin(locs{1}(1))*cos(locs{2}(1))*cos(delta_lon));
    course = 180/pi*theta;
    
    
end