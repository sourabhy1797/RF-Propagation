classdef Location
    %LOCATION Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        Name = '';
        Latitude;
        Longitude;
    end
    
    methods
        function obj = Location(lat,long)
            %LOCATION Construct an instance of this class
            %   Detailed explanation goes here
            obj.Longitude = long;
            obj.Latitude = lat;
        end
    end
end

