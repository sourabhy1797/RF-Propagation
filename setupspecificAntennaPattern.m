function setupspecificAntennaPattern(tx, angle, size)
    lambda = physconst("lightspeed") / tx.TransmitterFrequency;
    tx.Antenna  = phased.URA('Size',size,...
        'Lattice','Rectangular','ArrayNormal','x');
    % The multiplication factor for lambda units to meter conversion
    tx.Antenna .ElementSpacing = [0.5 0.5]*lambda/2;
  
    numRows = size(1);
    numCols = size(2);

    % Calculate Row taper
    rwind = ones(1,numRows).';
    % Calculate Column taper
    cwind = ones(1,numCols).';
    % Calculate taper
    taper = rwind*cwind.';
    tx.Antenna .Taper = taper.';
    % 
    % % % Create an sinc antenna element
    Elem = phased.SincAntennaElement;
    % Gaussian Antenna Element
    % Elem = phased.GaussianAntennaElement;
    Elem.FrequencyRange = [0 tx.TransmitterFrequency];
    Elem.Beamwidth = [10 10];
    tx.Antenna.Element = Elem;
    tx.AntennaAngle = angle;

end
