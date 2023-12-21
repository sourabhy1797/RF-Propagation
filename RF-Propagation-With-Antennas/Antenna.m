classdef Antenna
    %TRANSMITTER Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        Name
        Sim
        Type
        OpFreq
        Height = 0;
        Power = 10;
        Sensitivity
        MaxRange
        Tiltx = 0;
        Tilty = 0;
    end
    
    methods
        
        function obj = Antenna(Name, Type, OpFreq)
            %TRANSMITTER Construct an instance of this class
            %   Detailed explanation goes here
            obj.Name = Name;
            obj.Type = Type;
            obj.OpFreq = OpFreq;
            obj.Sim = design(Type, OpFreq);
        end

        function output = Generate_T_Site(this, name, location)
        output = txsite("Name",name , ...
                "Latitude", location.Latitude, ...
                "Longitude", location.Longitude, ...
                "Antenna", this.Sim, ...
                "AntennaHeight", this.Height, ...
                "TransmitterFrequency", this.OpFreq, ...
                "TransmitterPower", this.Power);
        end

        function output = Generate_R_Site(this, name, location)
        output = rxsite("Name",name , ...
                "Latitude", location.Latitude, ...
                "Longitude", location.Longitude, ...
                "Antenna", this.Sim, ...
                "ReceiverSensitivity", this.Sensitivity);
        end
    end
end

