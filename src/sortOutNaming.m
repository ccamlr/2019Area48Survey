function [stratum, transect] = sortOutNaming(t)
% Given a combined stratum and transect label, return the stratum name and
% transect label separately.

    if startsWith(t, "AP")
        stratum = "AP"; % Antarctic Pennisula
        transect = extractAfter(t,2);
    elseif startsWith(t, "ELE") % AMLR stratum
        stratum = "Elephant";
        transect = extractAfter(t,3);
    elseif startsWith(t, "WEST")  % AMLR stratum
        stratum = "West";
        transect = extractAfter(t,4);
    elseif startsWith(t, "SSI")
        stratum = "SSI";
        transect = extractAfter(t,3);
    elseif startsWith(t, "SSA") || startsWith(t, "SSB") || startsWith(t, "SSC")
        stratum = "ESS";
        transect = extractAfter(t,2);
    elseif startsWith(t, "SS") % Scotia Sea. Must come after SSI
        stratum = "SS";
        transect = extractAfter(t,2);
    elseif startsWith(t, "SG")
        stratum = "SG";
        transect = extractAfter(t,2);
    elseif startsWith(t, "SOF") % South Orkneys Fixed
        stratum = "SOF";
        transect = extractAfter(t,3);
    elseif startsWith(t, "SOC") % South Orkney Concentrated
        stratum = "SOC";
        transect = extractAfter(t,3);
    elseif startsWith(t, "SO") % South Okrney Islands. Must come after SOF and SOC
        stratum = "SOI";
        transect = extractAfter(t,2);
    elseif startsWith(t, "SA483")
        stratum = "SA483";
        transect = extractAfter(t,6);
    elseif startsWith(t, "Sand")
        stratum = "Sand";
        transect = extractAfter(t,4);
    elseif startsWith(t, "WCB")
        stratum = "WCB";
        transect = extractAfter(t,3);
    else
        stratum = "Unknown";
        transect = "Unknown";
    end
end