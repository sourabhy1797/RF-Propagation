classdef Region
    %REGION Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        NorthEastCorner
        SouthWestCorner
    end
    
    methods
        function obj = Region(NorthEastCorner, SouthWestCorner)
            %REGION Construct an instance of this class
            %   Detailed explanation goes here
            obj.NorthEastCorner = NorthEastCorner;
            obj.SouthWestCorner = SouthWestCorner;
        end
        
        function outputArg = North(obj)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            outputArg = obj.NorthEastCorner.Latitude;
        end

        function output = South(obj)
            output = obj.SouthWestCorner.Latitude;
        end

        function output = West(obj)
            output = obj.SouthWestCorner.Longitude;
        end

        function output = East(obj)
            output = obj.NorthEastCorner.Longitude;
        end
    end
end

