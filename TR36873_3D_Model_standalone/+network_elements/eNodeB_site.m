classdef eNodeB_site < handle
    % This object represents the eNodeB site with properties
    %
    % id        ... the eNodeB site id starting from 1
    % pos       ... position in [m]
    % altitude  ... altitude in [m]
    % site_type ... site type in case of a heterogeneus network (macro, pico, femto)
    % sectors   ... sector properties of the eNodeB site
    %
    % (c) Fjolla Ademaj, Martin Taranetz, ITC 2016
    properties
        id
        pos
        altitude
        site_type
        sectors
    end
    
end

