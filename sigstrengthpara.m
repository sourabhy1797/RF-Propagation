function ss = sigstrengthpara(rxs, txs, varargin)
%sigstrength   Received signal strength
%   SS = sigstrength(RX,TX) returns the signal strength at receiver site RX
%   due to transmitter site TX. 
%
%   SS = sigstrength(RX,TX,PROPMODEL) returns the signal strength with the
%   PropagationModel argument set to PROPMODEL. The default value
%   depends on the coordinate system in use by the input sites:
%
%      CoordinateSystem | Default PropagationModel value
%      --------------------------------------------------------------------
%      'geographic'     | 'longley-rice' when terrain is in use or else 
%                       | 'freespace' when terrain is not in use
%      --------------------------------------------------------------------
%      'cartesian'      | 'freespace' when Map is 'none' or else 
%                       | 'raytracing' when Map is an STL file or 
%                       | triangulation object or cartesian siteviewer
%
%   SS = sigstrength(___,Name,Value) returns the signal strength with
%   additional options specified by one or more Name-Value pairs.
%
%   The inputs RX and TX can be scalars or arrays. The output SS is an
%   array of size M-by-N, where M is the number of sites in TX and N is the
%   number of sites in RX.
%
%   sigstrength Name-Value pairs:
%
%   Type - Type of signal strength to compute, specified as 'power' (the
%      default value) or 'efield'. When Type is 'power', SS is expressed in
%      power units (dBm) of signal at the receiver input. When Type is
%      'efield', SS is expressed in electric field strength units (dBuV/m)
%      of signal wave incident on the antenna.
%
%   PropagationModel - Propagation model for path loss calculation,
%      specified as an object created with propagationModel or as one of
%      the following options:
%
%         'freespace'       - Free space propagation model
%         'rain'            - Propagation model for rain
%         'gas'             - Propagation model for gas
%         'fog'             - Propagation model for fog
%         'close-in'        - Close-In propagation model
%         'longley-rice'    - Longley-Rice propagation model
%         'tirem'           - Terrain Integrated Rough Earth Model (TIREM)
%         'raytracing'      - Ray tracing propagation model
%
%      Terrain propagation models, including the 'longley-rice' and 'tirem'
%      options, are only supported for sites with CoordinateSystem set to
%      'geographic'. The default value is described in the main help above.
%
%   Map - Map for surface geometry data. Valid and default values depend on
%      the coordinate system in use by the input sites:
%
%      Coordinate
%      system       | Valid Map values              | Default Map value
%      --------------------------------------------------------------------
%      'geographic' | siteviewer object or terrain  | current siteviewer
%                   | name ('none', 'gmted2010', or | object or else
%                   | the name of custom terrain    | 'gmted2010' if none
%                   | data previously added with    | is open
%                   | addCustomTerrain)             |
%      --------------------------------------------------------------------
%      'cartesian'  | 'none', triangulation object, | current siteviewer or
%                   | or the name of an STL file    | 'none'
%
%   % Example 1: Calculate received power and link margin at receiver
%
%   % Create transmitter site
%   tx = txsite('Name','Fenway Park', ...
%       'Latitude', 42.3467, ...
%       'Longitude', -71.0972);
%
%   % Create receiver site with sensitivity defined (in dBm)
%   rx = rxsite('Name','Bunker Hill Monument', ...
%       'Latitude', 42.3763, ...
%       'Longitude', -71.0611, ...
%       'ReceiverSensitivity', -90);
%
%   % Calculate received power and link margin, which is the difference
%   % between the receiver's sensitivity and the received power 
%   ss = sigstrength(rx,tx);
%   margin = abs(rx.ReceiverSensitivity - ss);
%
%   % Example 2: Compare received power of two sites in a conference room
%   % using ray tracing
%       
%   % Launch site viewer with a conference room 3D map
%   siteviewer("SceneModel", "conferenceroom.stl");
% 
%   % Define a transmitter site near the ceiling. Show it in site viewer.
%   tx = txsite("cartesian", ...
%       "AntennaPosition", [1.3; 1.3; 2.5], ...
%       "TransmitterFrequency", 2.8e9);
%   show(tx);
%
%   % Define two receiver sites: one above the table and the other below
%   % the table. Show them in site viewer. 
%   rxs = rxsite("cartesian", ...
%       "AntennaPosition", [0.3 0.3; 0.2 0.2; 1.0 0.6]);
%   show(rxs); 
% 
%   % Create a ray tracing propagation model for Cartesian coordinates 
%   % using the SBR method with up to 3 reflections. Set the surface 
%   % material to wood.
%   pm = propagationModel("raytracing", ...
%       "CoordinateSystem", "cartesian", ...
%       "Method","sbr", ...
%       "AngularSeparation", "low", ...
%       "MaxNumReflections", 3, ...
%       "SurfaceMaterial", "wood");
% 
%   % Calculate the received power at the two receiver sites from the
%   % transmitter site, using the ray tracing propagation model. The first
%   % receiver site has a much larger received power because it has a
%   % line-of-sight path to the transmitter. 
%   ss = sigstrength(rxs, tx, pm)
%   
%   See also link, raytrace, propagationModel

%   Copyright 2017-2022 The MathWorks, Inc.   

% Validate sites
validateattributes(rxs,{'rxsite'},{'nonempty'},'sigstrength','',1);
validateattributes(txs,{'txsite'},{'nonempty'},'sigstrength','',2);


% Allocate output matrix
numTx = numel(txs);
numRx = numel(rxs);
ss = zeros(numTx, numRx);

% Process optional inputs
p = inputParser;
if nargin > 2 && mod(numel(varargin),2)
    % Validator function is necessary for inputParser to allow string
    % option instead of treating it like parameter name
    p.addOptional('PropagationModel', [], @(x)ischar(x)||isstring(x)||isa(x,'rfprop.PropagationModel'));
else
    p.addParameter('PropagationModel', []);
end
p.addParameter('Type', 'power');
p.addParameter('ReceiverGain', []);
p.addParameter('Map', []);
p.addParameter('TransmitterAntennaSiteCoordinates', []);
p.addParameter('ReceiverAntennaSiteCoordinates', []);
p.addParameter('TerrainProfiles', []);
p.parse(varargin{:});

% Get usingCartesian from CoordinateSystem validation or from pre-validated 
% AntennaSiteCoordinates
if isempty(p.Results.TransmitterAntennaSiteCoordinates)
    usingCartesian = rfprop.internal.Validators.validateCoordinateSystem(rxs, txs);
else
    usingCartesian = strcmp(p.Results.TransmitterAntennaSiteCoordinates.CoordinateSystem,'cartesian');
end

% Validate and get parameters
if usingCartesian
    map = rfprop.internal.Validators.validateCartesianMap(p);
    pm = rfprop.internal.Validators.validateCartesianPropagationModel(p, map, 'sigstrength');
    if isprop(pm, 'CoordinateSystem')
        pm.CoordinateSystem = 'cartesian';
    end
else
    map = rfprop.internal.Validators.validateMapTerrainSource(p, 'sigstrength');
    rfprop.internal.Validators.validateMapCoordinateSystemMatch(map, 'geographic');
    pm = rfprop.internal.Validators.validateGeographicPropagationModel(p, map, 'sigstrength');
end
rfprop.internal.Validators.validateMaxNumReflections(pm, 'sigstrength');
isMultipathModel = pm.isMultipathModel;
typeIsEfield = validateType(p);
if typeIsEfield
    % Assuming far-field conditions and reciprocity:  
    %   radiated power density Wt (or St) = Pt*Gt*Lt*Fp^2/(4*pi*r^2);
    %       Pr = Wt*Gr*Lr*Fpol^2*lambda^2/(4*pi) [Friss Transmission Equation];
    %       and Wt = |E x conj(H)|avg = |Erms|^2 / eta. Fp is the
    %       propagation factor accounting for scenario complexity (=1 for
    %       homogenous medium e.g. antennas in free-space) and eta is the
    %       characteristic impedance of the propagation medium.
    % pathloss PL = (4*pi*r/(lambda*Fp))^2
    %   => Fp^2/(4*pi*r^2) = 4*pi/(PL*lambda^2)
    % => Erms = sqrt(eta*Wt) = sqrt(eta*Pt*Gt*4*pi/(PL*lambda^2))
    %   Note:  <ray>.PathLoss includes polarization mismatch (Fpol^2) so
    %   the electric field strength calculated herein is actually that
    %   which would be measured by the specified receive antenna (including
    %   its orientation).
    Z0 = rfprop.Constants.Z0; % propagation medium assumed to be free-space
    EfieldConst = 10*log10(Z0*4*pi) + 120; % dBuV/m offset constant
end

rxGain = validateReceiverGain(p, numTx);

% Get site antenna coordinates
txsCoords = rfprop.internal.Validators.validateAntennaSiteCoordinates(...
    p.Results.TransmitterAntennaSiteCoordinates, txs, map, 'sigstrength');
rxsCoords = rfprop.internal.Validators.validateAntennaSiteCoordinates(...
    p.Results.ReceiverAntennaSiteCoordinates, rxs, map, 'sigstrength');

% Calculate path loss depending on propagation model type
if isMultipathModel
    rays = raytrace(txs, rxs, pm, "Map", map);
    Lpls = rays;
    rxsAoas = rays;
else
    % Get path loss and propagation path info
    [Lpls_db, info] = pm.pathloss(rxs, txs, 'Map', map, ...
        'TransmitterAntennaSiteCoordinates', txsCoords, ...
        'ReceiverAntennaSiteCoordinates', rxsCoords, ...
        'ComputeAngleOfArrival', isempty(rxGain), ...
        'TerrainProfiles', p.Results.TerrainProfiles);
end

% Optimize pattern computation with re-usable directivity objects when
% transmitters share frequency and all tx antennas or all rx antennas are
% the same Phased Array object.
isSharedTxFrequency = all(arrayfun(@(x)x.TransmitterFrequency==txs(1).TransmitterFrequency, txs(:)));
sharedTxAntenna = txs(1).Antenna;
useSharedTxDirectivity = (numTx > 1) && ~isMultipathModel && isSharedTxFrequency && ...
    rfprop.AntennaSite.isPhasedAntenna(sharedTxAntenna) && ...
    all(arrayfun(@(x)x.Antenna==sharedTxAntenna, txs(:)));
if useSharedTxDirectivity
    txDirectivity = phased.internal.Directivity('Sensor',sharedTxAntenna);
end

sharedRxAntenna = rxs(1).Antenna;
useSharedRxDirectivity = (numRx > 1) && ~isMultipathModel && isSharedTxFrequency && ...
    rfprop.AntennaSite.isPhasedAntenna(sharedRxAntenna) && ...
    all(arrayfun(@(x)x.Antenna==sharedRxAntenna, rxs(:)));
if useSharedRxDirectivity
    rxDirectivity = phased.internal.Directivity('Sensor',sharedRxAntenna);
end

% Compute signal strength from each transmitter to all receivers
for txInd = 1:numTx
    tx = txs(txInd);

    % Compute transmitter constants
    fq = tx.TransmitterFrequency;
    Ptx = tx.TransmitterPower;
    if typeIsEfield
        lambda = rfprop.Constants.LightSpeed/fq;
        txEfieldConst = EfieldConst + 10*log10(Ptx) - 20*log10(lambda); % dBuV/m offset constant
    else
        Ptx_db = 10 * log10(1000*Ptx); % Convert W to dBm (db with reference to mW)
    end
    Ltxsys_db = tx.SystemLoss;

    % Get directivity or gain pattern
    useTxDirectivity = useSharedTxDirectivity || rfprop.AntennaSite.isPhasedAntenna(tx.Antenna);
    useTxGainPattern = rfprop.AntennaSite.isElectromagneticAntenna(tx.Antenna);
    if useTxDirectivity && ~useSharedTxDirectivity
        txDirectivity = phased.internal.Directivity('Sensor',tx.Antenna);
    elseif useTxGainPattern
        [Gtx,Gtxaz,Gtxel] = gainPattern(tx,fq);
    end

    % Get AoD for all propagation paths from tx. For multipath models e.g.
    % ray-trace, also need AoA for all tx-rx pairs and corresponding path
    % loss; path loss will be kept in dB if only 1 ray found per site pair,
    % otherwise path loss will be a unitless field strength.
    if isMultipathModel
        txAods = [];
        numRays = zeros(1,numRx);
        parfor rxInd = 1:numRx
            rxRays = rays{txInd,rxInd};
            if ~isempty(rxRays)
                rxNumRays = numel(rxRays);
                numRays(rxInd) = numel(rxRays);
                if rxNumRays == 1
                    txAods = [txAods rxRays.AngleOfDeparture]; %#ok<AGROW>
                    rxsAoas{txInd,rxInd} = rxRays.AngleOfArrival; % keep in dB
                    Lpls{txInd,rxInd} = rxRays.PathLoss;
                else
                    rxAoas_temp = zeros(2,rxNumRays);
                    Lpls_temp = zeros(rxNumRays,1);
                    for i = 1:rxNumRays
                        txAods = [txAods rxRays(i).AngleOfDeparture]; %#ok<AGROW>
                        rxAoas_temp(:,i) = rxRays(i).AngleOfArrival;
                        Lpls_temp(i) = exp(1i*rxRays(i).PhaseShift)/ ...
                            10^(rxRays(i).PathLoss/20);
                    end
                    rxsAoas{txInd,rxInd} = rxAoas_temp;
                    Lpls{txInd,rxInd} = Lpls_temp; % unitless field strength
                end
            end
            
        end
    else
        txAods = [info(txInd,:).AngleOfDeparture];
    end

    % Compute gains for all propagation paths from tx
    if ~isempty(txAods)
        txAz = txAods(1,:)';
        txEl = txAods(2,:)';
        if useTxDirectivity
            Gtxrxs_db = gain(tx,fq,txAz,txEl,txDirectivity);
        elseif useTxGainPattern
            Gtxrxs_db = gain(tx,fq,txAz,txEl,Gtx,Gtxaz,Gtxel);
        else
            Gtxrxs_db = gain(tx,fq,txAz,txEl);
        end
    end

    aodStartInd = 1;
    for rxInd = 1:numRx
        % Compute transmitter gain (including system losses) and path loss
        if isMultipathModel
            % Get path loss and check for empty, which corresponds to no
            % propagation paths found
            Lpl_ul = Lpls{txInd,rxInd};
            if isempty(Lpl_ul) % Infinite path loss since no propagation path
                ss(txInd,rxInd) = -inf;
                continue
            end

            % Get gains for this rx
            aodEndInd = aodStartInd + numRays(rxInd) - 1;
            Gtx_db = Gtxrxs_db(aodStartInd:aodEndInd) - Ltxsys_db;
            aodStartInd = aodEndInd + 1;
        else
            % Get path loss
            Lpl_db = Lpls_db(txInd, rxInd);

            Gtx_db = Gtxrxs_db(rxInd) - Ltxsys_db;
        end

        % Compute signal strength depending on if calculating E-Field
        % strength or received power. For E-Field strength, do not need to
        % know rx gain; for received power, do need to know rx gain.
        if typeIsEfield % dBuv/m
            % For multipath model e.g. ray-trace, use coherent phasor sum
            % of the individual path losses.
            if isMultipathModel
                % Perform phasor sum. If only 1 ray was found for this site
                % pair, the path loss will still be in dB, otherwise the
                % path loss will be unitless field strength.
                if numRays(rxInd) == 1 % Note: == 0 case is previously handled
                    ss(txInd,rxInd) = txEfieldConst + Gtx_db - Lpl_ul; % Lpl_ul is still in dB for only 1 ray
                else
                    E = sum(10.^(Gtx_db./20).*Lpl_ul);
                    ss(txInd,rxInd) = txEfieldConst + 20*log10(abs(E));
                end
            else
                ss(txInd,rxInd) = txEfieldConst + Gtx_db - Lpl_db; 
            end
        else % dBm
            % Compute receiver gain, including system losses
            if ~isempty(rxGain)
                Grx_db = repmat(rxGain(txInd), size(Gtx_db));
            else
                rx = rxs(rxInd);

                % Get directivity or gain pattern
                useRxDirectivity = useSharedRxDirectivity || rfprop.AntennaSite.isPhasedAntenna(rx.Antenna);
                useRxGainPattern = rfprop.AntennaSite.isElectromagneticAntenna(rx.Antenna);
                if useRxDirectivity && ~useSharedRxDirectivity
                    rxDirectivity = phased.internal.Directivity('Sensor',rx.Antenna);
                elseif useRxGainPattern
                    [Grx,Grxaz,Grxel] = gainPattern(rx,fq);
                end

                if isMultipathModel
                    rxAoas = rxsAoas{txInd,rxInd};
                else
                    infoStructs = info(txInd,rxInd);
                    rxAoas = [infoStructs.AngleOfArrival];
                end

                rxAz = rxAoas(1,:)';
                rxEl = rxAoas(2,:)';

                if useRxDirectivity
                    Grx_db = gain(rx,fq,rxAz,rxEl,rxDirectivity) - rx.SystemLoss;
                elseif useRxGainPattern
                    Grx_db = gain(rx,fq,rxAz,rxEl,Grx,Grxaz,Grxel) - rx.SystemLoss;
                else
                    Grx_db = gain(rx,fq,rxAz,rxEl) - rx.SystemLoss;
                end
            end

            % Compute signal strength in dBm, using link budget form of
            % Friis equation:
            %
            %  Received Power (dBm) = Transmitted Power (dBm) + Gains (dB) 
            %  - Losses (dB)
            %
            % Applied to tx/rx, this yields:
            %
            %  Prx_db = Ptx_db + Gtx_db + Grx_db - Lpl_db
            %
            % where:
            %  * Prx_db is received power in dBm at receiver input
            %  * Ptx_db is transmitter output power in dBm
            %  * Gtx_db is transmitter system gain in dBi (antenna gain - system loss)
            %  * Grx_db is receiver system gain in dBi (antenna gain - system loss)
            %  * Lpl_db is path loss (dB) as given by propagation model
            %
            % For multipath model e.g. ray-trace, use coherent phasor sum
            % of the individual path losses.
            if isMultipathModel
                % Perform phasor sum. If only 1 ray was found for this site
                % pair, the path loss will still be in dB, otherwise the
                % path loss will be unitless field strength.
                if numRays(rxInd) == 1 % 0 case is already taken care of
                    ss(txInd,rxInd) = Ptx_db + Gtx_db + Grx_db - Lpl_ul; % Lpl_ul is still in dB for only 1 ray
                else
                    E = sum(10.^(Gtx_db./20).*Lpl_ul.*10.^(Grx_db./20));
                    ss(txInd,rxInd) = Ptx_db + 20*log10(abs(E));
                end
            else
                ss(txInd,rxInd) = Ptx_db + Gtx_db + Grx_db - Lpl_db;
            end
        end
    end
end
end

function typeIsEfield = validateType(p)

try
    type = p.Results.Type;
    type = validatestring(type, {'efield','power'}, 'sigstrength','Type');
    typeIsEfield = strcmp(type,"efield");
catch e
    throwAsCaller(e);
end
end

function rxGain = validateReceiverGain(p, numTx)

try
    rxGain = p.Results.ReceiverGain;
    if ~isempty(rxGain)
        validateattributes(rxGain,{'numeric'}, ...
            {'real','finite','nonnan','nonsparse','nonempty'}, 'sigstrength', 'ReceiverGain');
        
        % Expand scalar gain to match length of tx
        if isscalar(rxGain)
            rxGain = repmat(rxGain,1,numTx);
        end
    end
catch e
    throwAsCaller(e);
end
end
