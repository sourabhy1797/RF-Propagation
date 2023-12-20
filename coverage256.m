function varargout = coverage256(txs, varargin)
%coverage   Display or compute coverage map
%   coverage(TX) displays the coverage map for transmitter site TX in the
%   current Site Viewer. Each colored contour of the map defines the area
%   where the corresponding signal strength is transmitted to a mobile
%   receiver. The map contours are generated using signal strength values
%   computed for receiver site locations on the map.
%
%   coverage(TX,PROPMODEL) displays the coverage map with the
%   PropagationModel argument set to PROPMODEL. The default propagation
%   model is 'longley-rice' when terrain is in use or else 'freespace' 
%   when terrain is not in use.
%
%   coverage(TX,RX) displays the coverage map with receiver arguments
%   ReceiverGain and ReceiverAntennaHeight set using receiver site RX.
%
%   coverage(TX,RX,PROPMODEL) displays the coverage map with receiver
%   arguments set using receiver site RX and the PropagationModel argument
%   set to PROPMODEL.
%
%   coverage(___,Name,Value) displays the coverage map with additional
%   options specified by one or more Name-Value pairs.
%
%   PD = coverage(TX,___) returns computed signal coverage data in
%   propagation data object PD. No plot is displayed and any graphical-only
%   Name-value pairs are ignored.
%
%   The input TX can be scalar or an array. If an array, the coverage map
%   displays the maximum signal strength received from any transmitter
%   site.
%
%   This function only supports antenna sites with CoordinateSystem set to
%   'geographic'.
%
%   coverage Name-Value pairs:
%
%   Type - Type of coverage map, specified as 'power' (the default value)
%      or 'efield'. When Type is 'power', SignalStrengths is expressed as
%      power units (dBm) of signal at the mobile receiver input. When Type
%      is 'efield', SignalStrengths is expressed in electric field strength
%      units (dBuV/m) of signal wave incident on the antenna.
%
%   SignalStrengths - Signal strengths to display in coverage map,
%      specified as a numeric vector. Each strength is displayed as a
%      different colored, filled contour on a map. Default value is -100
%      dBm if Type is 'power' and 40 dBuV/m if Type is 'efield'.
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
%      The default value is described in the main help above.
%
%   MaxRange - Maximum range of coverage map from each transmitter site,
%      specified as a positive numeric scalar in meters representing great
%      circle distance. MaxRange defines the region of interest on the map
%      to plot. The default value is automatically computed based on the 
%      propagation model type as shown:
%
%      Propagation model type   | MaxRange
%      ------------------------------------------------------------
%      Atmospheric or Empirical | Range of minimum value in SignalStrengths
%      Terrain                  | 30 km or distance to furthest building
%      Ray Tracing              | 500 m
%
%   Resolution - Resolution of receiver site locations used to compute
%      signal strength values, specified as 'auto' (the default value) or a
%      numeric scalar. The resolution defines the maximum distance between
%      locations. A Resolution of 'auto' computes the value based on
%      MaxRange. A numeric Resolution is expressed as a distance in meters.
%      Decreasing the Resolution value increases both the quality of the
%      coverage map and the time required to create it.
%
%   ReceiverGain - Mobile receiver's gain in dB, specified as a numeric
%      scalar. The ReceiverGain is used to compute received signal strength
%      when Type is 'power'. ReceiverGain includes the mobile receiver's
%      antenna gain and system loss. If receiver site argument RX is passed
%      to coverage, the default value is the maximum gain of the RX Antenna
%      with SystemLoss subtracted; otherwise the default value is 2.1.
%
%   ReceiverAntennaHeight - Mobile receiver's antenna height above ground 
%      elevation in meters, specified as a numeric scalar. If receiver site
%      argument RX is passed to coverage, the default value is the
%      AntennaHeight of RX; otherwise the default value is 1.
%
%   Colors - Colors of coverage map filled contours, specified as M-by-3
%      array of RGB (red, green, blue) triplets or array of strings or cell
%      array of character vectors. Colors are assigned element-wise to
%      values in SignalStrengths for coloring the corresponding filled
%      contour. Colors may not be used with Colormap or ColorLimits.
%
%   ColorLimits - Color limits for colormap, specified as a two-element
%      vector of the form [min max]. The color limits indicate the signal
%      level values that map to the first and last colors in the colormap.
%      Default value is [-120 -5] if Type is 'power' and [20 135] if Type
%      is 'efield'. ColorLimits may not be used with Colors.
%
%   Colormap - Colormap for coloring coverage map filled contours,
%      specified as a predefined colormap name or an M-by-3 array of RGB
%      (red, green, blue) triplets that define M individual colors. Default
%      value is 'jet'. Colormap may not be used with Colors.
%
%   ShowLegend - Show signal strength color legend on map, specified as
%      true or false. Default value is true.
%
%   Transparency - Transparency of coverage map, specified as a numeric
%      scalar in range 0 to 1, where 0 is completely transparent and 1 is
%      opaque. Default value is 0.4.
%
%   Map - Map for visualization and surface data, specified as a siteviewer
%      object or a terrain name. A terrain name may be specified if the
%      function is called with an output argument. Valid terrain names are
%      'none', 'gmted2010', or the name of custom terrain data previously
%      added with addCustomTerrain. The default value is the current
%      siteviewer. If no siteviewer is open, the default value is a new
%      siteviewer or else 'gmted2010' if the function is called with an
%      output argument.
%
%   % Examples
%   
%   % Example 1: Display coverage map for strong and weak signals
%   
%   % Define strong and weak signal strengths with corresponding colors
%   strongSignal = -70; % Unit: dBm
%   strongSignalColor = "green";
%   weakSignal = -90; % Unit: dBm
%   weakSignalColor = "cyan";
%
%   % Create transmitter site and display coverage map
%   tx = txsite("Name","MathWorks", ...
%       "Latitude", 42.3001, ...
%       "Longitude", -71.3503, ...
%       "AntennaHeight", 50);
%   coverage(tx, ...
%      "SignalStrengths",[strongSignal,weakSignal], ...
%      "Colors", [strongSignalColor,weakSignalColor])
%
%   % Update coverage map with new signal strengths. The coverage map is
%   % updated instead of added since the same txsite object is used.
%   strongSignal = -65; % Unit: dBm
%   weakSignal = -80; % Unit: dBm
%   coverage(tx, ...
%      "SignalStrengths",[strongSignal,weakSignal], ...
%      "Colors", [strongSignalColor,weakSignalColor])
%
%   % Example 2: Display coverage map with buildings
%   
%   % Launch Site Viewer with buildings in Chicago
%   viewer = siteviewer("Basemap","streets-light",...
%      "Buildings","chicago.osm");
%
%   % Create transmitter site on a building
%   tx = txsite("Latitude",41.8800, ...
%      "Longitude",-87.6295, ...
%      "TransmitterFrequency",2.5e9);
%   show(tx)
%
%   % Define a ray tracing propagation model using the
%   % shooting-bouncing-rays (SBR) method up to three reflections
%   pm = propagationModel("raytracing", ...
%       "Method", "sbr", ...
%       "MaxNumReflections", 3);
%
%   % Display coverage map using the ray tracing propagation model, which
%   % computes total received power as the sum of power along all
%   % propagation paths. Limit coverage area to 200 m from txsite.
%   coverage(tx,pm,...
%      "SignalStrengths",-90:-5,...
%      "Resolution",5,...
%      "MaxRange",200,...
%      "Transparency",0.6)
%
%   % Example 3: Display combined coverage map for multiple transmitters
%
%   % Define names and locations of sites around Boston
%   names = ["Fenway Park","Faneuil Hall","Bunker Hill Monument"];
%   lats = [42.3467,42.3598,42.3763];
%   lons = [-71.0972,-71.0545,-71.0611];
%
%   % Create transmitter site array
%   txs = txsite("Name", names,...
%      "Latitude",lats,...
%      "Longitude",lons, ...
%      "TransmitterFrequency",2.5e9);
%
%   % Display combined coverage map for multiple transmitter sites, using
%   % Close-In propagation model. The map shows the maximum signal strength
%   % received from any transmitter site.
%   coverage(txs,"close-in", ...
%      "SignalStrengths",-100:5:-60)
%
% See also <a href="matlab:help rxsite.sigstrength">sigstrength</a>, <a href="matlab:help txsite.sinr">sinr</a>, <a href="matlab:help rxsite.link">link</a>, propagationModel

%   Copyright 2017-2022 The MathWorks, Inc.

% Validate number of output arguments
nargoutchk(0,1)

% Validate site
validateattributes(txs,{'txsite'},{'nonempty'},'coverage','',1);

% Validate sites are all geographic
rfprop.internal.Validators.validateGeographicSites(txs, 'coverage')

% Process optional rxsite input
rx = [];
numTx = numel(txs);
args = varargin;
defaultRxGain = 2.1;
defaultRxHeight = 1;
if numel(args) > 0
    secondInput = args{1};
    if isa(secondInput,'rxsite')
        rx = secondInput;
        validateattributes(rx,{'rxsite'},{'scalar'},'coverage','',2);
        args = args(2:end); % Remove from input list
        
        % Set default receiver parameters from input. Generate default rx
        % gain for each tx, since they may have different frequencies.
        defaultRxGain = zeros(1,numTx);
        for k = 1:numTx
            defaultRxGain(k) = gain(rx,txs(k).TransmitterFrequency) - rx.SystemLoss;
        end
        defaultRxHeight = rx.AntennaHeight;
    end
end

% Add parameters
p = inputParser;
if mod(numel(args),2) % Odd number of inputs
    % Validator function is necessary for inputParser to allow string
    % option instead of treating it like parameter name
    p.addOptional('PropagationModel', [], @(x)ischar(x)||isstring(x)||isa(x,'rfprop.PropagationModel'));
else
    p.addParameter('PropagationModel', []);
end
p.addParameter('Type', 'power');
p.addParameter('SignalStrengths', []);
p.addParameter('Resolution', 'auto');
p.addParameter('ReceiverGain', defaultRxGain);
p.addParameter('ReceiverAntennaHeight', defaultRxHeight);
p.addParameter('Animation', '');
p.addParameter('MaxRange', 'auto');
p.addParameter('Colormap', 'jet');
p.addParameter('Colors', [])
p.addParameter('ColorLimits', []);
p.addParameter('Transparency', 0.4);
p.addParameter('ShowLegend', true);
p.addParameter('ReceiverLocationsLayout', []);
p.addParameter('MaxImageResolutionFactor', 5);
p.addParameter('RadialResolutionFactor', 2);
p.addParameter('Map', []);
p.parse(args{:});

% Get Site Viewer and validate web graphics
outputRequested = nargout > 0;
if ~outputRequested
    viewer = rfprop.internal.Validators.validateMap(p, 'coverage');
    isViewerInitiallyVisible = viewer.Visible;
else
    isViewerInitiallyVisible = false;
end

% Get coverage map type and corresponding values 
allowedTypes = {'efield','power'};
[type, defaultStrengths, defaultColorLimits, legendTitle] = ...
    rfprop.internal.Validators.validateType(p, allowedTypes, 'coverage');

% Validate parameters
map = rfprop.internal.Validators.validateMapTerrainSource(p, 'coverage');
terrainSource = rfprop.internal.Validators.validateTerrainSource(map, 'coverage');
animation = rfprop.internal.Validators.validateGraphicsControls(p, isViewerInitiallyVisible, 'coverage');
strengths = validateSignalStrengths(p, defaultStrengths);
colors = rfprop.internal.Validators.validateColors(p, 'coverage');
clim = rfprop.internal.Validators.validateColorLimits(p, defaultColorLimits, 'coverage');
cmap = rfprop.internal.Validators.validateColorMap(p, 'coverage');
transparency = rfprop.internal.Validators.validateTransparency(p, 'coverage');
showLegend = rfprop.internal.Validators.validateShowLegend(p, 'coverage');
pm = rfprop.internal.Validators.validateGeographicPropagationModel(p, map, 'coverage');
rfprop.internal.Validators.validateMaxNumReflections(pm, 'coverage');
rfprop.internal.Validators.validateMaxNumDiffractions(pm, 'coverage'); % Disable 2 order diffractions
rxGain = rfprop.internal.Validators.validateReceiverGain(p, 'coverage');
rxAntennaHeight = rfprop.internal.Validators.validateReceiverAntennaHeight(p, 'coverage');
maxImageResFactor = p.Results.MaxImageResolutionFactor;
radialResFactor = p.Results.RadialResolutionFactor;

% Get site coordinates
txsCoords = rfprop.internal.AntennaSiteCoordinates.createFromAntennaSites(txs, map);
txslatlon = txsCoords.LatitudeLongitude;

% Validate dependent parameters
if isscalar(rxGain) % Expand scalar gain to match length of tx
    rxGain = repmat(rxGain,1,numTx);
end
maxrange = validateMaxRange(p, type, pm, txs, txslatlon, strengths, rxGain, map);
[res, isAutoRes] = rfprop.internal.Validators.validateResolution(p, maxrange, 'coverage');
datarange = rfprop.internal.Validators.validateDataRange(txslatlon, maxrange, res, ~strcmp(terrainSource,'none'));
rxLocationsLayout = rfprop.internal.Validators.validateReceiverLocationsLayout(p, pm, txslatlon, 'coverage');

if outputRequested
    % Do not show progress dialog
    generatingMapMsg = '';
    computingDataMsg = '';
else
    % Get color info
    colorsSpecified = ~ismember('Colors', p.UsingDefaults);
    colormapSpecified = ~all(ismember({'Colormap','ColorLimits'}, p.UsingDefaults));
    useColors = colorsSpecified || (numel(strengths) < 2 && ~colormapSpecified);
    if useColors
        if isempty(colors) && ismember(type,{'power','efield'})
            colors = [0 1 0]; % Default coverage map color (green)
        end
        colorData = struct('Levels',strengths,'Colors',colors);
    else
        colorData = struct('Colormap',cmap,'ColorLimits',clim);
    end
    
    % Clear old contour color info for each site. This is because a site can
    % never conflict with itself (since a site will replot it's own contour if
    % coverage is called on a site multiple times)
    % Also mark the contour graphics to be removed
    graphicsToRemove = {};
    for k = 1:numTx
        oldContourID = viewer.getGraphic(txs(k).UID, 'contour');
        if (~isempty(oldContourID))
            viewer.removeColorData(oldContourID);
            graphicsToRemove = [graphicsToRemove; oldContourID]; %#ok<AGROW>
        end
    end
    if ~isempty(viewer.LegendID)
        graphicsToRemove = [graphicsToRemove; viewer.LegendID];
    end
    
    % Generate a new ID so old IDs don't conflict with new images
    idNum = viewer.getId(1);
    contourID = ['contour' num2str(idNum{1})];
    
    % Validate image for color conflicts. This must be done before calculations
    % occur to avoid wasting the user's time when color conflict occurs.
    oldColorGraphics = viewer.ColorGraphics;
    viewer.checkForGraphicsConflict(type, contourID, colorData);
    resetColorGraphics = rfprop.internal.onExit(@()setColorGraphics(viewer,oldColorGraphics));
    
    % Show txsites. If Site Viewer is already open, keep current camera view.
    if isViewerInitiallyVisible && ~viewer.Visible
        return % Abort if Site Viewer has been closed
    end
    if isViewerInitiallyVisible
        showAnimation = 'none';
    else
        showAnimation = 'zoom';
    end
    show(txs,'Map',viewer,'Animation',showAnimation, ...
        'AntennaSiteCoordinates', txsCoords);
    
    % Show progress dialog. Use cleanup object to close dialog if a forced
    % exit occurs (error or Ctrl-C). Use three-stage dialog if using ray
    % tracing and single-stage dialog if not a terrain model. If using a
    % terrain model, dialog is launched in radialReceiverLocationsLayoutData.
    resetProgressDialog = rfprop.internal.onExit(@()hideProgressDialog(viewer));
    showWaitbar = pm.isMultipathModel;
    generatingMapMsg = message('shared_channel:rfprop:ProgressDialogGeneratingCoverageMap').getString;
    computingDataMsg = message('shared_channel:rfprop:ProgressDialogComputingCoverage').getString;
    if ~pm.requiresTerrain || strcmp(rxLocationsLayout,'grid')
        if showWaitbar
            msg = message('shared_channel:rfprop:ProgressDialogPreparingMapData').getString;
        else
            msg = generatingMapMsg;
        end
        viewer.showProgressDialog('Message', msg, ...
            'Indeterminate', true, ...
            'Cancelable', false);
    end
end

% Generate location grid containing data range from each transmitter site
if isa(map,'siteviewer')
    maxImageSize = map.MaxImageSize;
else
    maxImageSize = rfprop.Constants.DefaultMaxImageSize;
end
[latNorth, latSouth, lonEast, lonWest, animation] = ...
    rfprop.internal.MapUtils.geobounds(txslatlon, datarange, animation);
[gridlats, gridlons, res] = rfprop.internal.MapUtils.geogrid(...
    latNorth, latSouth, lonEast, lonWest, res, isAutoRes, maxrange, maxImageSize, 'coverage');
gridSize = size(gridlats);

% Compute and validate coverage map image size
imageSize = rfprop.internal.Validators.validateImageSize(...
    gridSize, maxImageResFactor, maxImageSize, res, 'coverage');

gridlats = gridlats((length(gridlats) / 2) - 128 : (length(gridlats) / 2) + 128 , (length(gridlats) / 2) - 128 : (length(gridlats) / 2) + 128);
gridlons = gridlons((length(gridlons) / 2) - 128 : (length(gridlons) / 2) + 128 , (length(gridlons) / 2) - 128 : (length(gridlons) / 2) + 128);

% Trim grid locations to those which are within data range
[latitude, longitude] = rfprop.internal.MapUtils.georange(...
    txs, gridlats(:), gridlons(:), datarange, terrainSource);

% Trim grid locations to those which are not inside buildings
if pm.isMultipathModel && isa(map,'siteviewer')
    isWithinBldg = isWithinAnyBuilding(map.BuildingsArray, latitude, longitude);
    latitude = latitude(~isWithinBldg);
    longitude = longitude(~isWithinBldg);
end

if ~outputRequested && siteViewerWasClosed(viewer)
    return % Abort if Site Viewer has been closed
end
if isempty(latitude)
    warning(message('shared_channel:rfprop:NoCoverageMapArea'))
    if outputRequested
        varargout = {[]};
    else
        viewer.removeColorData(viewer.SiteGraphics.(txs(k).UID).contour);
        viewer.remove(graphicsToRemove);
    end
    return
end

% Use try/catch in case cancel error thrown from waitbar dialog
try
    if strcmp(rxLocationsLayout,'grid')
        % Define rxsites at grid locations within data range
        rxs = rxsite(...
            'Name', 'internal.coveragesite', ... % Specify to avoid default site naming
            'Latitude', latitude, ...
            'Longitude', longitude, ...
            'AntennaHeight', rxAntennaHeight);
        
        % Launch waitbar phase of dialog
        if ~outputRequested && showWaitbar
            viewer.showProgressDialog('Message', computingDataMsg, ...
                'Value', 0, ...
                'ShowPercentage', true, ...
                'Indeterminate', false, ...
                'Cancelable', true);
        end
        
        % Compute signal strength at each rxsite in the grid.
        data = sigstrength(rxs, txs, pm, ...
            'Type', type, ...
            'ReceiverGain', rxGain, ...
            'Map', map, ...
            'TransmitterAntennaSiteCoordinates', txsCoords);
        
    else
        % Compute data and corresponding locations using radial layout
        [latitude, longitude, ~, data] = radialReceiverLocationsLayoutData(txs, txsCoords, type, ...
            generatingMapMsg, computingDataMsg, pm, map, res, radialResFactor, datarange, ...
            lonWest, lonEast, latSouth, latNorth, rxGain, rxAntennaHeight);
    end
catch e
    if strcmp(e.identifier, "shared_channel:rfprop:ProgressDialogCancelled")
        return
    else
        rethrow(e)
    end
end

% Launch final phase of status dialog 
if ~outputRequested
    viewer.showProgressDialog('Message', generatingMapMsg, ...
        'Indeterminate', true, ...
        'Cancelable', false);
end

% Merge data of rx sites for use in single image, where signal strength at
% an rx site is the max strength due to any tx.
if numTx > 1
    data = max(data, [], 1);
end

% Remove NaN-data
isnn = ~isnan(data);
data = data(isnn);
latitude = latitude(isnn);
longitude = longitude(isnn);

if ~outputRequested 
    % Show rxsite
    if siteViewerWasClosed(viewer)
        return
    end
    if isa(rx,'rxsite')
        show(rx,'Map',viewer,'Animation','none','EnableWindowLaunch',false);
    end
end

% Create propagation data container
if strcmpi(type,'power')
    dataVariableName = 'Power';
else
    dataVariableName = 'Efield'; 
end
pd = propagationData(latitude, longitude, ...
    'Name', message('shared_channel:rfprop:CoveragePropagationDataTitle').getString, ...
    dataVariableName, data(:));

if outputRequested
    % Return propagationData with contour properties set except for ID,
    % which is not set since there is no plot to associate with
    pd = pd.setContourProperties(txslatlon(:,1),txslatlon(:,2),maxrange,imageSize);
    varargout = {pd};
    return
end

% It is now safe to remove old contours associated with sites
viewer.remove(graphicsToRemove);

% Return early if no coverage map data meets minimum signal strength,
% including if no rays were found (all data = -Inf)
minstrength = min(strengths,[],"all");
[mindata, maxdata] = bounds(data,"all");
if maxdata < minstrength
    if maxdata == -Inf % No rays found
        warning(message("shared_channel:rfprop:NoCoverageMapAreaHasFiniteSignalStrength"))
    else
        warning(message("shared_channel:rfprop:NoCoverageMapAreaMeetsMinSignalStrength", ...
            num2str(minstrength),num2str(mindata),num2str(maxdata)))
    end
    
    return
end

% Create contour map
pd = pd.setContourProperties(txslatlon(:,1),txslatlon(:,2),maxrange,imageSize,contourID);
contourArgs = {'Type', type, ...
    'Map', viewer, ...
    'Levels', strengths, ...
    'Transparency', transparency, ...
    'ShowLegend', showLegend, ...
    'ImageSize', imageSize, ...
    'LegendTitle', legendTitle, ...
    'ValidateColorConflicts', false};
if useColors
    contourArgs = [contourArgs, ...
        'Colors', colors];
else
    contourArgs = [contourArgs, ...
        'ColorLimits', clim, ...
        'Colormap', cmap];
end
contour(pd, contourArgs{:});

% Set the IDs of the appropriate Site Graphics
for i = 1:numel(txs)
    viewer.SiteGraphics.(txs(i).UID).contour = contourID;
    if (showLegend)
        viewer.SiteGraphics.(txs(i).UID).legend = ['legend' contourID];
    end
end

% Hide progress dialog
viewer.hideProgressDialog;

% Cancel forced exit cleanup, since normal execution has completed
cancel(resetProgressDialog);
if ~outputRequested
    cancel(resetColorGraphics);
end
end

function wasClosed = siteViewerWasClosed(viewer)

wasClosed = viewer.LaunchWebWindow && ~viewer.Visible;
end

function strengths = validateSignalStrengths(p, defaultStrengths)

try
    if ismember('SignalStrengths', p.UsingDefaults)
        strengths = defaultStrengths;
    else
        strengths = p.Results.SignalStrengths;
        validateattributes(strengths, {'numeric'}, ...
            {'real','finite','nonnan','nonsparse','vector','nonempty'}, 'coverage', 'SignalStrengths');
        if ~iscolumn(strengths)
            strengths = strengths(:);
        end
        strengths = double(strengths);
    end
catch e
    throwAsCaller(e);
end
end

function maxrange = validateMaxRange(p, type, pm, txs, txslatlon, strengths, rxGain, map)

try
    maxrange = p.Results.MaxRange;
    numSites = numel(txs);
    if ischar(maxrange) || isstring(maxrange)
        validatestring(maxrange, {'auto'}, 'coverage', 'MaxRange');
        
        % Get max range for each tx
        if pm.requiresTerrain || pm.isMultipathModel
            % Use default value for terrain prop model
            maxrange = rfprop.internal.Validators.defaultMaxRange(txslatlon, pm, map);
        else
            maxrange = zeros(1,numSites);
            ss = min(strengths);
            
            % Turn off warning on calling range since check is made below
            warnState = warning('off','shared_channel:rfprop:RangeGreaterThanMax');
            warnCleanup = onCleanup(@()warning(warnState));
            for k = 1:numSites
                maxrange(k) = range(txs(k), ss, rxGain(k), type, pm);
            end
            
            % Get propagation limit
            if isa(map,'siteviewer') 
                useTerrain = map.UseTerrain;
            else
                useTerrain = ~strcmp(map,'none');
            end
            if useTerrain
                maxRangeLimit = rfprop.Constants.MaxPropagationDistanceUsingTerrain;
                warnID = 'shared_channel:rfprop:MaxRangeAssignedTerrainLimit';
            else
                maxRangeLimit = rfprop.Constants.MaxPropagationDistance;
                warnID = 'shared_channel:rfprop:MaxRangeAssignedPropagationLimit';
            end
            
            % Saturate and throw warning if range is at limit. Use limit
            % equality instead of strict greater than check since "range"
            % above may have already saturated to the limit.
            exceedsMax = (maxrange >= maxRangeLimit);
            if any(exceedsMax)
                maxrange(exceedsMax) = maxRangeLimit;
                warning(message(warnID,round(maxRangeLimit/1000)));
            end
        end
    else
        maxrange = rfprop.internal.Validators.validateNumericMaxRange(maxrange, pm, numSites, map, 'coverage');
    end
catch e
    throwAsCaller(e);
end
end

function isWithinBldg = isWithinAnyBuilding(bldgs, lats, lons)

% Check if locations are within footprint of any buildings
isWithinBldg = false(numel(lons),1);
for bldgInd = 1:numel(bldgs)
    [inp, onp] = bldgs(bldgInd).Footprint.isinterior(lats,lons);
    isWithinBldg = isWithinBldg | inp | onp;
end
end

function setColorGraphics(viewer,imageGraphics)

viewer.ColorGraphics = imageGraphics;
end
