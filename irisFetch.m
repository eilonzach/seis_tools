classdef irisFetch
   %IRISFETCH allows seamless access to data stored within the  IRIS-DMC
   %
   %Each method retrieves data from the IRIS-DMC, returned in array of
   %matlab structures.  To see the available methods, type:
   %
   %  methods(irisFetch)
   %
   %
   %OVERVIEW OF METHODS:
   %
   %  IRISFETCH.Traces() retrieves sac-equivelent waveform with channel
   %  metadata from the IRIS-DMC
   %
   %    Example:
   %      tr = IRISFETCH.Traces('IU','ANMO','10','BHZ',...
   %           '2010-02-27 06:30:00','2010-02-27 10:30:00')
   %
   %
   %  IRISFETCH.Stations() retrieves station metadata from the IRIS-DMC.
   %  This data may be retrieved at a variety of detail levels. From
   %  broadest to most specific, these are NETWORK, STATION, CHANNEL, and
   %  RESPONSE.
   %
   %    Example:
   %      s = IRISFETCH.Stations('channel','*','ANMO','10','BH?')
   %      % the above returns a tree-like structure, now make it more manageable
   %      s = IRISFETCH.flattenToChannel(s); % convert a 1xN array of channel epochs.
   %
   %
   %  IRISFETCH.Events() retrieves event metadata from the IRIS-DMC
   %    Example:
   %      ev = IRISFETCH.Events('MinimumMagnitude',6.0,...
   %          'minimumLatitude',45,'maximumLatitude', 60,...
   %          'minimumLongitude', -150,'maximumLongitude', -90)
   %
   %
   %
   %REQUIREMENTS:
   %
   %  IRISFETCH requires the latest version of the Web Services Library,
   %  which is a java .jar file available from:
   %
   %  http://www.iris.edu/manuals/javawslibrary/#download
   %
   %  This jar file must be added to your MATLAB path, which may be done
   %  in a variety of ways.  One common way is to include a javaaddpath
   %  statement in the startup.m file.  For more details, consult MATLAB's
   %  documentation for 'Bringing Java Classes and Methods into MATLAB
   %  Workspace'.
   %
   %
   %FOR FURTHER GUIDANCE:
   %
   %  A discussion about using MATLAB to access data from the IRIS-DMC can
   %  be found at:
   %
   %  http://www.iris.edu/manuals/javawslibrary/matlab/
   %
   %
   %
   %see also IRISFETCH.TRACES IRISFETCH.STATIONS IRISFETCH.EVENTS
   %IRISFETCH.FLATTENTOCHANNEL IRISFETCH.FLATTENTOSTATION
   %IRISFETCH.TESTCOMPATIBILITY IRISFETCH.VERSION JAVAADDPATH
   
   
   % Celso Reyes, Rich Karstens
   % IRIS-DMC
   % February 2012
   
   % 2012 Nov 7, r1.3.6
   % Fixed an occasional rounding error at the ms level.
   % Duplicated PreferredMagnitude and PreferredOrigin information into the main level of
   % the structure returned by Events.  This should make dealing with the data much
   % easier.
   %
   % 2012 Sept 26 r1.3.5
   % Refractored to make irisFetch easier to read
   % added some more error catching codes, in effort to make program
   % debugging easier.  Added support for RESP webservice (library version 1.5.XXX)
   %
   % 2012 July 19 r1.3.4
   % Fixed error for Stations with responses that have different types of
   % responses.  Depending upon the version, matlab throws out two
   % different (but similar) error ids.
   %
   % 2012 June 18 r1.3.2
   % spelling fix in parse, and initialized
   %
   % 2012 June 14 r1.3.1
   % Fixed problem where Traces.sacpz.units was not converted from java
   % strings
   %
   % 2012 June 8, r1.3.0
   % Changed Traces routine to match IRIS-WS v1.5's new conventions.
   %  - This change should be transparent to matlab users (except new
   %    version of library must be downloaded)
   %  - verbosity set through routine, not as parameter
   %  - fetch sacpz is required
   % added sacpz subfield: units.
   % fixed bug in Event where radius didn't convert from double to
   % java.lang.Double, resulting in error.
   % added the App name to Traces.
   % added flattenToStation routine
   %
   % Modified the Traces routine: Added authorization ability and ability to get
   % poles and zeros responses.  Reworked the parameter list to be more
   % flexible. Completely overhauled the java->matlab parsing for accuracy and speed.
   %
   % 2012 Mar 8, Fixed problem where empty traces created a "data not
   % assigned" style message. Removed java objects from the returned
   % structures... these have already been parsed into fields.
   %
   % 2012 Feb 29, Fixed date parsing issues, allowing accuracy to
   % millisecond level. Also fixed problem where filter numerators not
   % showing up properly.
   %
   % 2012 Feb 27, 1.1.1 fixed problem where channel epochs were missing from
   % flattened stations, and added ability to translate into both Long and
   % Integer java types.
   %
   % 2012 Feb 13, added SampleRate to trace structure, improved error
   % catching.
   %
   % 2012 Feb 9, minor documentation update
   
   
   properties (Constant = true)
      VERSION = '1.3.5';
      MS_IN_DAY = 86400000;
      DATE_FORMATTER = 'yyyy-mm-dd HH:MM:SS.FFF';
      BASE_DATENUM = 719529; % accounts for matlab's 0000-Jan-1 start date vs java's 1970-Jan-1 start
      MIN_JAR_VERSION = '1.5';
      SURROGATE_JAR = 'http://www.iris.edu/manuals/javawslibrary/matlab/IRIS-WS-1.5-matlab.jar';
      
      VALID_QUALITIES = {'D','R','Q','M','B'};
      DEFAULT_QUALITY = 'B';
      
   end %properties
   
   methods(Static)
      function v = version()
         % return the version number of irisFetch
         v = irisFetch.VERSION;
      end
      
      
      
      function ts = Traces(network, station, location, channel, startDateStr, endDateStr, varargin )
         %irisFetch.Traces accesses waveform with associated channel metadata
         %
         %USAGE
         % tr = irisFetch.Traces(network, station, location,
         % channel, startDate, endDate) will use channel and date
         % criteria to retrieve one or more seismic traces, which are
         % stored as structures containing typical SAC-equivalent
         % metadata.
         %
         % startDate and endDate must be formatted thusly:
         %      'YYYY-MM-DD hh:mm:ss' or 'YYYY-MM-DD hh:mm:ss.sss'
         %
         % network, station, location, and channel all accept '*' and
         % '?' wildcards, as well as comma separated lists.
         %
         % tr = irisFetch.Traces(..., quality)
         % allows you to specify a quality factor, such as 'B' or 'M'.
         % If not specified, quality defaults to 'B'
         %
         % tr = irisFetch.Traces(..., 'includePZ')
         % will retrieve poles and zeroes as well.  These are the same
         % poles/zeros as found from http://www.iris.edu/ws/sacpz/
         %
         % tr = irisFetch.Traces(..., 'verbose')
         % provides additional debugging information
         %
         % tr = irisFetch.Traces(..., usernameAndPassword )
         % allows authorized users to access restricted data.
         % usernameAndPassword must be a cell containing the username and
         % password.
         %    Sample:
         %      unamepwd = {'myname@iris.edu', 'mypassword'}
         %
         %
         %
         %EXAMPLES:
         %
         %    Example 1:
         %      % get 3 components for a 4-minute time period,
         %      % using the '?' wildcard and also retrieving response info.
         %      ts = irisFetch.Traces('IU','ANMO','10','BH?',...
         %           '2010-02-27 06:30:00','2010-02-27 10:30:00','includePZ')
         %
         %    Example 2:
         %      % get z-channels for a comma-separated list of stations
         %      % while specifying a quality of 'B'
         %      ts = irisFetch.Traces('IU','ANMO,ANTO,YSS','00','BHZ',...
         %           '2010-02-27 06:30:00','2010-02-27 10:30:00', 'B')
         %
         %    Example 3:
         %      % get z-data for all BHZ stations that belong to the IU
         %      % network and have a location code of '00'
         %      ts = irisFetch.Traces('IU','*','00','BHZ',...
         %           '2010-02-27 06:30:00','2010-02-27 10:30:00')
         %
         %ABOUT THE RETURNED TRACE
         %  The returned trace(s) will be an array of structures, with
         %  each structure containing fields full of information. For
         %  example, if I retrieve N traces, I will have the following
         %  structure:
         %
         %  1xN struct array with fields:
         %     network
         %     station
         %     location
         %     channel
         %     quality
         %     latitude
         %     longitude
         %     elevation
         %     depth
         %     azimuth
         %     dip
         %     sensitivity
         %     sensitivityFrequency
         %     instrument
         %     sensitivityUnits
         %     data
         %     sampleCount
         %     sampleRate
         %     startTime
         %     endTime
         %
         %   COMMON MANIPULATIONS
         %     To access the date as text, use datestr().  For trace(s)
         %     ts the full date with milliseconds can be seen using:
         %     datestr(ts, ,'YYYY-MM-DD hh:mm:ss.FFF')
         %
         %     To scale the data, divide each trace's data by its
         %     sensitivity.  The resulting values are in
         %     sensitivityUnits.
         %
         % SEE ALSO datestr
         
         % these variables are shared among all the nested functions
         getsacpz = false;
         verbosity = false;
         authorize = false;
         quality = irisFetch.DEFAULT_QUALITY;
         username = '';
         userpwd = '';
         tracedata = []; % initializing for scoping issues
         
         extractAdditionalArguments(varargin);
         conformifyDates();
         conformifyLocation();
         connectToIrisLibrary();
         setAppName();
         setVerbosity();
         getTheTraces();
         return
         
         % ---------------------------------------------------------------
         % END TRACES: MAIN
         % ===============================================================
         
         
         % ---------------------------------------------------------------
         % TRACES: NESTED FUNCTIONS
         % v---v---v---v---v---v---v---v---v---v---v---v---v---v---v---v
         
         function extractAdditionalArguments(argList)
            % extracts getsacpz, verbosity, authorize, quality, username, and userpwd
            % Parameters are handled "intelligently" so that [paramname, paramval] pairs
            % aren't necessry
            
            for n=1:numel(argList)
               thisParameter = argList{n};
               switch class(thisParameter)
                  case 'cell'
                     setCredentials(thisParameter)
                     authorize = true;
                  case 'char'
                     switch upper(thisParameter)
                        case irisFetch.VALID_QUALITIES
                           quality = thisParameter;
                        case {'INCLUDEPZ'}
                           getsacpz = true;
                        case {'VERBOSE'}
                           verbosity = true;
                        otherwise
                           error('IRISFETCH:Trace:unrecognizedParameter',...
                              'The text you included as an optional parameter did not parse to either a qualitytype (D,R,Q,M,B) or ''INCLUDEPZ'' or ''VERBOSE''');
                     end
                  case 'logical'
                     verbosity = thisParameter; % old usage, may be deprecated in the future.
                  otherwise
                     error('IRISFETCH:Trace:unrecognizedParameter','The optional parameter wasn''t recognized. %s', class(thisParameter));
               end
            end
            
            % debug display
            % disp({'spz:',getsacpz,'vb:',verbosity,'auth:',authorize,'qual:',quality,'un&pw:',username,userpwd}); % DEBUG
            
            function setCredentials(param)
               % parameter is the username/pwd combo
               assert(numel(param)==2 && ischar(param{1}) && ischar(param{2}),...
                  ['A cell passed as an optional parameter is assumed',...
                  ' to contain credentials. eg. {''myname'',''mypassword''}.']);
               username = param{1};
               userpwd = param{2};
            end % setCredentials
            
            
         end % extractAdditionalArguments
         
         function conformifyDates()
            startDateStr = irisFetch.makeDateStr(startDateStr);
            endDateStr = irisFetch.makeDateStr(endDateStr);
         end %conformifyDates
         
         function conformifyLocation()
            location = strrep(location,' ','-');
         end %conformifyLocaiton
         
         function connectToIrisLibrary()
            try
               tracedata = edu.iris.dmc.ws.extensions.fetch.TraceData();
            catch er
               switch er.identifier
                  case {'MATLAB:undefinedVarOrClass', 'MATLAB:subscripting:undefinedClass'}
                     warning('IRISFETCH:NoIrisWSJarInstalled',MessageDownloadLibrary());
                     if ~irisFetch.connectTo_IRIS_WS_jar('silent')
                        error('IRISFETCH:Traces:UnableToInstallIrisWSJar',...
                           'irisFetch was unable to recover, please download and add the latest IRIS-WS-JAR to your javaclasspath');
                     end
                     disp('irisFetch.connectTo_IRIS_WS_jar() has was able to connect you to the appropriate java library. Continuing...');
                     tracedata = edu.iris.dmc.ws.extensions.fetch.TraceData();
                  otherwise
                     rethrow(er)
               end
            end
            
         end %connectToIrisLibrary
         
         function setAppName()
            tracedata.setAppName(['MATLAB:irisFetch/' irisFetch.version()]);
         end %setAppName
         
         function setVerbosity()
            try   % only library 1.5 and greater will successfully do this
               tracedata.setVerbosity(verbosity);
            catch er
               % if the error is due to a bad library version, then recommend
               % updating the library.
               
               %otherwise
               rethrow(er);
            end
         end %setVerbosity
         
         
         function getTheTraces()
            traces = []; %will be structure of traces.
            try
               fetchTracesBasedOnAuthorization();
            catch je
               switch je.identifier
                  case {'MATLAB:undefinedVarOrClass', 'MATLAB:subscripting:undefinedClass'}
                     attemptToRecoverFromMissingLibrary();
                     
                  case 'MATLAB:Java:GenericException'
                     if containsURLNotFoundException(je)
                        throwTraceUrlNotFoundException(je);
                     end
                     rethrow(je)
                     
                  otherwise
                     fprintf('Exception occured in IRIS Web Services Library: %s\n', je.message);
                     rethrow(je)
                     
               end
            end
            
            ts = irisFetch.convertTraces(traces);
            clear traces
            
            
            function val = containsURLNotFoundException(je)
               % we got a 404 from somewhere. (based on ice.net)
               val = any(strfind(je.message,'URLNotFoundException'));
            end
            
            function throwTraceUrlNotFoundException(je)
               error('IRISFETCH:Trace:URLNotFoundException','Trace found no requested data. Instead, it ran in to the following error:\n%s', je.message);
            end
            
            function fetchTracesBasedOnAuthorization()
               if authorize
                  fetchTracesWithAuthorization();
               else
                  fetchTracesNormally();
               end
            end
            
            function fetchTracesWithAuthorization()
               traces = tracedata.fetchTraces(network, station, location, channel, ...
                  startDateStr, endDateStr, quality, getsacpz, username, userpwd);
            end
            
            function fetchTracesNormally()
               traces = tracedata.fetchTraces(network, station, location, channel, startDateStr, endDateStr, quality, getsacpz);
            end
            
            function attemptToRecoverFromMissingLibrary()
               warning('IRISFETCH:NoIrisWSJarInstalled',MessageDownloadLibrary());
               
               if irisFetch.connectTo_IRIS_WS_jar('silent')
                  disp('irisFetch.connectTo_IRIS_WS_jar() has was able to connect you to the appropriate java library. Continuing...');
                  fetchTracesBasedOnAuthorization();
               else
                  error('IRISFETCH:Traces:UnableToInstallIrisWSJar',...
                     'irisFetch was unable to recover, please download and add the latest IRIS-WS-JAR to your javaclasspath');
               end
            end
            
         end %function getTheTraces
         
         function msg = MessageDownloadLibrary()
            msg=['The Web Services library does not appear to be in the javaclasspath.\n',...
               'Please download the latest version from \n',...
               'http://www.iris.edu/manuals/javawslibrary/#download\n ',...
               'and then add it to your classpath. \n'];
         end %downloadLatestVersionMessage
         
         
      end % Traces
      
      
      function [networkStructure, urlParams] = Stations(detailLevel, network, station, location, channel, varargin)
         %irisFetch.Stations retrieves station metadata from IRIS-DMC
         %
         %USAGE:
         %  s = irisFetch.Stations(detail, network, station, location,
         %  Channel),  where detail is one of "NETWORK", "STATION",
         %  "CHANNEL", or "RESPONSE".  These five parameters are
         %  required for all queries, but may be wildcarded by using []
         %  for their values.
         %
         %  Network, station, location, and channel parameters are
         %  passed directly to the java library, so lists (separated by
         %  commmas) and wildcards (? and *) are accepted.
         %
         %  [s, myParams] = irisFetch.Stations( ... ) also returns the
         %  URL parameters that were used to make the query.
         %
         %  s = irisFetch.Stations( ... , param1, value1 [, ...]])
         %  allows any number of parameter-value pairs to be included in
         %  the selection criteria.  All of the StationCriteria 'set'
         %  methods are supported as parameter pairs, as are a couple
         %  special parameter shortcuts.
         %
         % To determine which parameters can be set,
         %    1) type the following:
         %         methods(edu.iris.dmc.ws.criteria.StationCriteria)
         %    2) look at the list:
         %         All methods starting with "set" are accessible via
         %         the parameter list.
         %
         %   Example:
         %     n = irisFetch.Stations(....,'endbefore',now, 'maxlongitude',-100)
         %
         %     which would invoke the following setters:
         %       stationCriteria.setEndBefore( now )
         %       stationCriteria.setMaxLongitude(-100)
         %
         %
         %
         %  Usable parameters are listed below.  For detailed
         %  descriptions of their effect and use, consult the station
         %  webservice webpage, available at:
         %
         %  http://www.iris.edu/ws/station/
         %
         %PARAMETER LIST (as of 9/27/2012)
         %  'MinimumLatitude', 'MaximumLatitude', 'MinimumLongitude',
         %  'MaximumLongitude', 'Latitude', 'Longitude',
         %  'MinimumRadius','MaximumRadius', 'StartAfter', 'EndAfter',
         %  'StartBefore', 'EndBefore', 'StartTime', 'EndTime',
         %  'UpdatedAfter'
         %
         %CONVENIENCE PARAMETERS
         %   'boxcoordinates'    : [minLat, minLon, maxLat, maxLon]
         %                           % use NaN as a wildcard
         %   'radialcoordinates' : 1x4 double :
         %                           [Lat, Lon, MaxRadius, MinRadius]
         %                           % MinRadius is optional
         %
         %
         %ADDITIONAL DETAILS: USING CELLS
         %  Network, station, location, and channel may be strings or
         %  cell arrays. Each element of the cell array is added to the
         %  search criteria as "or".  That is:
         %
         %      irisFetch.Stations(detail,'AB,BC,CD',...)
         %
         %  is equivalent to
         %
         %      irisFetch.Stations(detail,{'AB','BC','CD'},...)
         %
         %  However, the former makes one call to addNetwork() while
         %  the latter makes three calls to addNetwork(). The net effect
         %  may be the same, but the execution is different.
         %
         %
         %WORKING WITH THE RESULTS
         %  The results are returned in a structure tree, with the same
         %  hierarchy found in the StationXML.  To make this easier to
         %  work with within matlab, you can use the
         %  irisFetch.flattenToChannel routine to shuffle the results
         %  into a 1xN array of Channels.  This only works if the detail
         %  level was "Channel" or "Response".
         %
         %SEE ALSO IRISFETCH.FLATTENTOCHANNEL
         
         %END OF HELP
         
         %-------------------------------------------------------------
         % An alternate base URL can be specified by providing the
         % parameter pair (...,'BASEURL',alternateURL)
         %
         % This is mostly useful in testing, and is likely not relevent
         % for users
         %-------------------------------------------------------------
         
         % import edu.iris.dmc.*
         % import edu.iris.dmc.ws.station.model.*
         
         outputLevel = '';
         service = []; % will be a java service object
         criteria = []; % will be a java criteria object
         j_networks = []; % will be the returned java networks
         
         % =============================================================
         % STATIONS: MAIN
         % v---v---v---v---v---v---v---v---v---v---v---v---v---v---v---v
         verifyArguments(nargin);
         setOutputLevel();
         connectToStationService();
         setCriteria();
         fetchTheStations();
         convertStationsToMatlabStructs();
         returnTheUrlParams();
         return
         
         % -------------------------------------------------------------
         % END STATIONS: MAIN
         % =============================================================
         
         
         % -------------------------------------------------------------
         % STATIONS: NESTED FUNCTIONS
         % v---v---v---v---v---v---v---v---v---v---v---v---v---v---v---v
         function verifyArguments(nArgs)
            if nArgs==1 && strcmpi(detailLevel,'help')
               disp('HELP request recognized, but not implemented');
               return
            elseif nArgs < 5
               error('not enough arguments.%d',nArgs);
            end
         end %verifyArguments
         
         function setOutputLevel()
            try
               outputLevel = edu.iris.dmc.ws.criteria.OutputLevel.(upper(detailLevel));
            catch je
               switch je.identifier
                  case 'MATLAB:undefinedVarOrClass'
                     error('IRISFETCH:NoIrisWSJarInstalled',...
                        'The necessary IRIS-WS java library was not recognized or found. Please ensure it is on your javaclasspath');
                  case 'MATLAB:subscripting:classHasNoPropertyOrMethod'
                     error('IRISFETCH:invalidOutputLevel',...
                        'The selected outputLevel [''%s''] was not recognized.',...
                        upper(detailLevel));
                  otherwise
                     rethrow(je);
               end
            end
         end % setOutputLevel
         
         function connectToStationService()            
            serviceManager = edu.iris.dmc.ws.service.ServiceUtil.getInstance();            
            serviceManager.setAppName(['MATLAB:irisFetch/' irisFetch.version()])            
            setBaseUrlFromParameterList();
            removeParameter('BASEURL');
            return
            
            % - - - - - - - - - - - - - - -
            function setBaseUrlFromParameterList()
               baseUrl = getParameter('BASEURL');
               if ~isempty(baseUrl)                  
                  service = serviceManager.getStationService(baseUrl);
               else
                  service = serviceManager.getStationService();                  
               end % setBaseUrlFromParameterList                 
            end
                       
         end %connectToStationService()
         
         function removeParameter(s)
            [~, idx] = getParameter(s);
            varargin(idx * 2 -1 : idx* 2) = [];
         end
         
         function [p, i] = getParameter(s)
            i = find(strcmpi(parameterNames(),s),1,'first');
            p = parameterValues();
            p = p(i);
         end
         
         function pn = parameterNames()
            pn = varargin(1:2:end);
         end
         
         function pv = parameterValues()
            pv = varargin(2:2:end);
         end
         
         
         function setCriteria()
            criteria = edu.iris.dmc.ws.criteria.StationCriteria;
            
            %----------------------------------------------------------
            % Deal with the Station/Network/Channel/Location parameters
            % These are treated separately, as they're "add" & not "set"
            % Each may handle multiple strings (as a cell array)
            %----------------------------------------------------------
            
            criteria = irisFetch.addCriteria(criteria, network, 'addNetwork');
            criteria = irisFetch.addCriteria(criteria, station, 'addStation');
            criteria = irisFetch.addCriteria(criteria, location,'addLocation');
            criteria = irisFetch.addCriteria(criteria, channel, 'addChannel');            
            criteria = irisFetch.setCriteria(criteria, varargin);            
         end %setCriteria
         
         function fetchTheStations()
            try
               j_networks = service.fetch(criteria, outputLevel);
            catch je
               if strfind(je.message,'ServiceNotSupportedException')
                  error('IRISFETCH:ServiceNotSupportedByLibrary',...
                     'The IRIS-WS java library version doesn''t support the requested station service version');
               else
                  rethrow(je)
               end
            end
         end %fetchTheStations
         
         function convertStationsToMatlabStructs()
            networkStructure = irisFetch.parse(j_networks);
         end
         
         function returnTheUrlParams()
            if nargout == 2
               urlParams = criteria.toUrlParams;
            end
         end
         
      end %Stations
      
      %%
      function [events, urlParams] = Events(varargin)
         %irisFetch.Events retrieves event data from the IRIS-DMC
         %
         %USAGE:
         %  ev = irisFetch.Events(param1, value1 [, ...]) retrieves
         %  event data from the IRIS-DMC database, returning it as a
         %  matlab structure.  An arbitrary number of parameter-value
         %  pairs may be specified in order to narrow down the search
         %  results.
         %
         %  [ev, myParams] = irisFetch.Events( ... ) also returns the
         %  URL parameters that were used to make the query.
         %
         %  Usable parameters are listed below.  For detailed
         %  descriptions of their effect and use, consult the webservice
         %  webpage for events, available at:
         %
         %  http://www.iris.edu/ws/event/
         %
         %PARAMETER LIST (as of 9/27/2012)
         %  'EventId'
         %  'FetchLimit'
         %  'MinLatitude'
         %  'MaxLatitude'
         %  'MinLongitude'
         %  'MaxLongitude'
         %  'MinimumDepth'
         %  'MaximumDepth'
         %  'Latitude'
         %  'Longitude'
         %  'MinRadius'
         %  'MaxRadius'
         %  'StartTime'
         %  'EndTime'
         %  'UpdatedAfter'
         %  'MagnitudeType'
         %  'MinimumMagnitude'
         %  'MaximumMagnitude'
         %  'catalogContains'
         %  'contributorContains'
         %  'includeArrivals'
         %  'includeallMagnitudes'
         %  'preferredOnly'
         %
         %
         %CONVENIENCE PARAMETERS
         %   'boxcoordinates'    : [minLat, minLon, maxLat, maxLon]   % use NaN as a wildcard
         %   'radialcoordinates' : [Lat, Lon, MaxRadius, MinRadius]   % MinRadius is optional
         %
         
         %END OF HELP
         
         %-------------------------------------------------------------
         % An alternate base URL can be specified by providing the
         % parameter pair (...,'BASEURL',alternateURL)
         %
         % This is mostly useful in testing, and is likely not relevent
         % for users
         %-------------------------------------------------------------
         
         import edu.iris.dmc.*
         
         serviceManager = ws.service.ServiceUtil.getInstance();
         serviceManager.setAppName(['MATLAB:irisFetch/' irisFetch.version()]);
         
         
         indexOffsetOfBASEURL=find(strcmpi(varargin(1:2:end),'BASEURL'),1,'first') * 2;
         
         try
            baseURL = varargin{indexOffsetOfBASEURL};
         catch
            % don't do anything
         end
         
         if exist('baseURL','var')
            varargin(indexOffsetOfBASEURL-1:indexOffsetOfBASEURL) = [];
            service = serviceManager.getEventService(baseURL);
         else
            service = serviceManager.getEventService();
         end
         
         criteria = ws.criteria.EventCriteria;
         criteria = irisFetch.setCriteria(criteria, varargin);
         if nargout == 2
            urlParams = criteria.toUrlParams;
         end
         disp('fetching from IRIS-DMC')
         j_events = service.fetch(criteria);
         fprintf('\n\n%d events found *************\n\n',j_events.size);
         disp('parsing into MATLAB structures')
         %tic;events = irisFetch.parseCollection(j_events);toc
         %disp(events)
         for n=size(j_events):-1:1
            %tic;
            %fprintf('parsing event %d :   ', n);
            thisEvent = irisFetch.parse(j_events.get(n-1));
            if numel(thisEvent.PreferredMagnitude)
                          thisEvent.PrimaryMagnitudeType=thisEvent.PreferredMagnitude.Type;
            thisEvent.PrimaryMagnitudeValue=thisEvent.PreferredMagnitude.Value;
            else
               thisEvent.PrimaryMagnitudeType='';
               thisEvent.PrimaryMagnitudeValue=nan;
            end
            if numel(thisEvent.PreferredOrigin)
               thisEvent.PrimaryLatitude=thisEvent.PreferredOrigin.Latitude;
               thisEvent.PrimaryLongitude=thisEvent.PreferredOrigin.Longitude;
               thisEvent.PrimaryDepth=thisEvent.PreferredOrigin.Depth;
               thisEvent.PrimaryTime=thisEvent.PreferredOrigin.Time;
            else
               thisEvent.PrimaryLatitude=nan;
               thisEvent.PrimaryLongitude=nan;
               thisEvent.PrimaryDepth=nan;
               thisEvent.PrimaryTime='0000-00-00 00:00:00.000';
            end
               
            events(n) = thisEvent;
            %disp(toc)
         end
         % events = irisFetch.parse(j_events);
      end
      
      
      %%
      function success = connectTo_IRIS_WS_jar(isSilent)
         %irisFetch.connectTo_IRIS_WS_jar tries to set up the jar for
         %use within MATLAB for this session
         %
         %USAGE:
         %  success = connectTo_IRIS_WS_jar('silent')
         %
         %  This routine searches the javaclasspath for the IRIS-WS jar
         %  file. If it does not exist, then it will try to access the
         %  latest jar over the internet.
         
         isSilent = exist('isSilent','var') && strcmpi(isSilent,'silent');
         
         success = false;
         %Check for required jar file for winston
         try
            % store the java class path for inspection
            jcp = javaclasspath('-all');
            
         catch er
            disp('Java is not enabled on this machine.  The Web Services Library will not work.');
            return
         end
         
         RequiredFiles = {'IRIS-WS-'};
         
         introuble = false;
         
         for FN = RequiredFiles
            if isempty(strfind([jcp{:}],FN{1}))
               if ~isSilent
                  disp(['Missing ' FN{1}]);
               end
               introuble = true;
            end
         end
         
         if introuble
            if ~isSilent
               disp('please add the IRIS-WS-latest.jar file to your javaclasspath');
               disp('ex.  javaaddpath(''/usr/local/somewhere/IRIS-WS.jar'');');
            end
            
            surrogate_jar = irisFetch.SURROGATE_JAR;
            
            [~,success] = urlread(surrogate_jar);%can we read the .jar? if not don't bother to add it.
            if success
               javaaddpath(surrogate_jar);
            else
               warning('irisFetch:noDefaultJar',...
                  'Unable to access the default jar.  Please download and add the latest IRIS-WS-JAR to your javaclasspath.');
            end
         end;
         success = true;
      end
      
      
      function channelList = flattenToChannel(networkTree)
         %irisFetch.flattenToChannel flattens the structure returned by irisFetch.Stations
         %
         %
         %USAGE
         %  flatStruct = irisFetch.flattenToChannel(networkTree)
         %
         %This takes the hierarchy returned by irisFetch.Stations, and
         %returns a 1xN array of channels  (channel epochs, technically).
         %
         % networkTree is a nested structure of the following format
         %   network.Station.Epoch.Channel.Epoch.[etc]
         %
         % flatStruct is an array containing ALL channel epochs, along
         % with unique identifying information from the parent levels,
         % such as networkTree.code, networkTree.station.code, etc.
         %
         %WORKING with the channelList
         %
         % Example: grabbing all BHZ channels from network IU
         %   myChannels = channelList({strcmp(channelList.NetworkCode},'IU') & ...
         %                   strcmp(channelList.ChannelCode, 'BHZ')
         %
         % Example: Find out how many IU networked stations were
         %          retrieved
         %   sum(strcmp({channelList.NetworkCode},'IU'))
         
         
         % moving from:  network -> station -> station epoch ->
         % channel -> channel epoch -> etc.
         
         % moving to: flat channels
         
         % hard-coded.
         %first, loop through and get rid of excess fields
         
         ALPHABETIZE = false; %leave fields in alphabetical order, or reorder according to common-sense
         
         if ~isa(networkTree,'struct')
            error('Cannot Flatten a non-structure');
         end
         
         % shared itterators
         eachStation = [];
         eachNetwork = [];
         
         % shared lists
         stationCodes = [];
         stationSites = []; 
         channelCodes = [];
         locationCodes = [];
         
         %Translate from a tree structure to a flat channel list   
         for eachNetwork = 1 : numel(networkTree)
            moveEverythingInTreeToTopLevel();
         end %eachNetwork
         channelList = deal([networkTree.Channels]);
         
         % clean up the Channel List
         moveSensitivityUpToMainLevel();
         moveSensorUpToMainLevel();
         
         if ~ALPHABETIZE
            reorderForCoherency();
         end
         
         return
         
         % ----------------------------------------------------------------
         % Flatten to Channel : Nested functions
         % ----------------------------------------------------------------
         
         function moveEverythingInTreeToTopLevel()            
            grabStationCodesForThisNetworkTree();
            for eachStation = 1 : numel(networkTree(eachNetwork).Stations)
               storeCurrentStationSiteList()
               organizeAtStationEpochLevel();
               moveChannelsFromEpochsToStationsLevel();
            end %eachStation
            stripEpochsFieldFromStations();
            moveChannelsUpToNetworkLevel();
         end
         
         
         function storeCurrentStationSiteList()
            stationSites = {networkTree(eachNetwork).Stations(eachStation).Epochs.Site};
         end
         
         function moveChannelsFromEpochsToStationsLevel()
               networkTree(eachNetwork).Stations(eachStation).Channels = ...
                  deal([networkTree(eachNetwork).Stations(eachStation).Epochs.Channels]);
         end
         
         function  moveSensitivityUpToMainLevel()
            % Bring Sensitivity and Sensor up, since it is a 1x1 struct
            moveFieldUpToMainLevelThenDelete('Sensitivity');
         end
         
         function  moveSensorUpToMainLevel()
            moveFieldUpToMainLevelThenDelete('Sensor');
         end
         
         function  moveFieldUpToMainLevelThenDelete(fieldToMigrate)
            for n=1:numel(channelList)
               if isstruct(channelList(n).(fieldToMigrate))
                  sensor=channelList(n).(fieldToMigrate);
                  fn = fieldnames(sensor);
                  for m=1:numel(fn);
                     channelList(n).(fn{m}) = sensor.(fn{m});
                  end
               end
            end
            channelList = rmfield(channelList,{(fieldToMigrate)});
         end
         
         function grabStationCodesForThisNetworkTree()
            stationCodes = {networkTree(eachNetwork).Stations.Code};
         end
         
         function stripEpochsFieldFromStations()
            networkTree(eachNetwork).Stations = rmfield(networkTree(eachNetwork).Stations,'Epochs');
         end
         
         function moveChannelsUpToNetworkLevel()
            networkTree(eachNetwork).Channels = deal([networkTree(eachNetwork).Stations.Channels]);
         end
         
         function organizeAtStationEpochLevel()
            for eachStationEpoch = 1 : numel(networkTree(eachNetwork).Stations(eachStation).Epochs)
               if ~isstruct(networkTree(eachNetwork).Stations(eachStation).Epochs(eachStationEpoch).Channels)
                  continue; % nothing to do.
               end
               storeChannelCodes();
               storeLocationCodes();
               organizeAtChannelLevel()
               networkTree(eachNetwork).Stations(eachStation).Epochs(eachStationEpoch).Channels = [networkTree(eachNetwork).Stations(eachStation).Epochs(eachStationEpoch).Channels.Epochs];
               % now, the structure is
               % net -> sta -> epoch -> chan
            end %eachStationEpoch
            
            function storeChannelCodes()
               channelCodes = {networkTree(eachNetwork).Stations(eachStation).Epochs(eachStationEpoch).Channels.Code};
            end
            
            function storeLocationCodes()
               locationCodes = {networkTree(eachNetwork).Stations(eachStation).Epochs(eachStationEpoch).Channels.Location};
            end
            
            function organizeAtChannelLevel()
               for eachChannel = 1 : numel(networkTree(eachNetwork).Stations(eachStation).Epochs(eachStationEpoch).Channels)
                  theseEpochs = networkTree(eachNetwork).Stations(eachStation).Epochs(eachStationEpoch).Channels(eachChannel).Epochs;
                  
                  [theseEpochs.NetworkCode] = deal(networkTree(eachNetwork).Code);
                  [theseEpochs.NetworkDescription] = deal(networkTree(eachNetwork).Description);
                  
                  [theseEpochs.StationCode] = deal(stationCodes{eachStation});
                  [theseEpochs.ChannelCode] = deal(channelCodes{eachChannel});
                  [theseEpochs.LocationCode]= deal(locationCodes{eachChannel});
                  [theseEpochs.Site] = deal(stationSites{eachStationEpoch});
                  networkTree(eachNetwork).Stations(eachStation).Epochs(eachStationEpoch).Channels(eachChannel).Epochs = theseEpochs;
               end %eachChannel
               
               
            end
         end
         

         
         function reorderForCoherency()
            % now, reorder to make it visually coherent.
            descriptorstuff={'NetworkCode';'StationCode';'LocationCode';'ChannelCode';'NetworkDescription';'Site'};
            positionalstuff={'Latitude';'Longitude';'Elevation';'Depth';'Azimuth';'Dip'};
            otherstuff={'SampleRate';'StartDate';'EndDate'};
            fieldsattop=[descriptorstuff; positionalstuff; otherstuff];
            
            fn = fieldnames(channelList);
            
            fieldsattop = fieldsattop(ismember(fieldsattop,fn)); %ensure fields exist
            
            for n=1:numel(fieldsattop);
               fn(strcmp(fn,fieldsattop(n))) = [];
            end
            neworder = [fieldsattop; fn];
            channelList = orderfields(channelList, neworder);
         end
         
      end
      
            
      function flatStruct = flattenToStation(networkTree)
         %irisFetch.flattenToStation flattens the structure returned by irisFetch.Stations
         %
         %
         %USAGE
         %  flatStruct = irisFetch.flattenToStation(networkTree)
         %
         %This takes the hierarchy returned by irisFetch.Stations, and
         %returns a 1xN array of stations (station epochs, technically).
         %
         % networkTree is a nested structure of the following format
         %   network.Station.Epoch.[etc]
         %
         % flatStruct is an array containing ALL station epochs, along
         % with unique identifying information from the parent levels,
         % such as networkTree.code, networkTree.station.code, etc.
         
         if isempty(networkTree)
            flatStruct = networkTree;
            return
         end
         
         if isempty([networkTree.Stations])
            flatStruct = networkTree;
            return
         end
         
         alphabetize = false; %leave fields in alphabetical order, or reorder according to common-sense
         
         migrateNetworkDetailsToStationLevel();         
         migrateStationDetailsToEpochLevel();         
         removeChannelsIfEmpty();
         
         if ~alphabetize
            reorderForCoherency();
         end
         
         return
         
         % ----------------------------------------------------------------
         % Flatten to Station : Nested functions
         % ----------------------------------------------------------------
         
         function migrateNetworkDetailsToStationLevel()            
            for netIdx = 1: numel(networkTree)
               % add the Network code to the Stations
               [networkTree(netIdx).Stations.NetworkCode] = deal(networkTree(netIdx).Code);
               [networkTree(netIdx).Stations.NetworkDescription] = deal(networkTree(netIdx).Description);
            end
            flatStruct = [networkTree.Stations];
         end
         
         function migrateStationDetailsToEpochLevel()
            
            for staIdx = 1 : numel(flatStruct);
               [flatStruct(staIdx).Epochs.StationCode] = deal(flatStruct(staIdx).Code);
               [flatStruct(staIdx).Epochs.NetworkCode] = deal(flatStruct(staIdx).NetworkCode);
               [flatStruct(staIdx).Epochs.NetworkDescription] = deal(flatStruct(staIdx).NetworkDescription);
            end
            flatStruct = [flatStruct.Epochs];
         end
         
         function removeChannelsIfEmpty()            
            if isempty([flatStruct.Channels]);
               flatStruct = rmfield(flatStruct,'Channels');
            end
         end
         
         function reorderForCoherency()
            % now, reorder to make it visually coherent.
            descriptorstuff={'NetworkCode';'StationCode';'LocationCode';'ChannelCode';'NetworkDescription';'Site'};
            positionalstuff={'Latitude';'Longitude';'Elevation';'Depth';'Azimuth';'Dip'};
            otherstuff={'SampleRate';'StartDate';'EndDate'};
            fieldsattop=[descriptorstuff; positionalstuff; otherstuff];
            
            fn = fieldnames(flatStruct);
            
            fieldsattop = fieldsattop(ismember(fieldsattop,fn)); %ensure fields exist
            
            for n=1:numel(fieldsattop);
               fn(strcmp(fn,fieldsattop(n))) = [];
            end
            neworder = [fieldsattop; fn];
            flatStruct = orderfields(flatStruct, neworder);
         end
         
      end
      
      function [respstructures, urlparams] = Resp(network, station, location, channel, starttime, endtime)
         % retrieve the RESP information into a character string.
         % net, sta, loc, and cha are all required.
         % channels and locations may be wildcarded using either ? or *
         % starttime and endtime options may be ignored by using [] instead of a time.
         
         criteria = edu.iris.dmc.ws.criteria.RespCriteria();
         criteria.setNetwork(network);
         criteria.setStation(station);
         criteria.setLocation(location);
         criteria.setChannel(channel);
         if ~isempty(starttime)
            criteria.setStartTime(irisFetch.mdate2jdate(starttime));
         end
         if ~isempty(endtime)
            criteria.setEndTime(irisFetch.mdate2jdate(endtime));
         end
         urlparams = criteria.toUrlParams().toCharArray()';
         
         serviceManager = edu.iris.dmc.ws.service.ServiceUtil.getInstance();
         baseUrl = 'http://www.iris.edu/ws/resp/';
         serviceManager.setAppName(['MATLAB:irisFetch/' irisFetch.version()]);
         service = serviceManager.getRespService(baseUrl);
         %respstructures= service.fetch(criteria);
         %respstructures = char(respstructures);
         respstructures= service.fetch(criteria).toCharArray()';
         
         
      end
      
      function [js, je] = testResp(starttime, endtime)
         n=now;
         testThis('IU','ANMO','00','BHZ',[],[]);
         testThis('IU','ANMO','00','*',[],[]);
         testThis('IU','ANMO','*','BHZ',[],[]);
         testThis('IU','ANMO','00','BHZ',[],now-1);
         testThis('IU','ANMO','00','BHZ',n-600,[]);
         testThis('IU','ANMO','00','BHZ',n,n-1);
         testThis('IU','ANMO','?0','BHZ',n,n-1);
         testThis('IU','ANMO','00','B?Z',n,n-1);
         testThis('IU','ANMO','00','BHZ','9/20/2012','9/20/2012 03:00:00');
         %testThis('IU','*','00','BHZ',[],[]); %should fail
         %testThis('*','ANMO','00','BHZ',[],[]); % should fail
         if exist('starttime','var') && ~isempty(starttime)
            js = showTimeInDetail(starttime);
         end
         disp(' ');
         if exist('endtime','var') && ~isempty(endtime)
            je = showTimeInDetail(endtime);
         end
         
         function testThis(varargin)
            try
               [r, url] = irisFetch.Resp(varargin{:});
               whos r
               disp(['url: ', url]);
               
               parampairs={'network',varargin{1},...
                  'station',varargin{2},...
                  'location',varargin{3},...
                  'channel',varargin{4}}
               
               if ~isempty(varargin{5})
                  st = datestr(varargin{5},31);
                  st(11)='T';
                  parampairs = [parampairs, {'starttime',st}];
               end
               
               if ~isempty(varargin{6})
                  ed = datestr(varargin{6},31);
                  ed(11)='T';
                  parampairs = [parampairs, {'endtime',ed}];
               end
               
               [s,code]=urlread('http://www.iris.edu/ws/resp/query','get', parampairs);
               
               assert(strcmp(r,s));
            catch myerror
               warning('RESPTEST:failure',myerror.identifier);
            end
         end
         
         function javadateTime = showTimeInDetail(t)
            dv=datevec(t);
            dv(6) = ceil(dv(6) * 1000) / 1000;
            t = datenum(dv);
            s = dv(6);
            % t must be either a date number or a string.
            javadateTime = irisFetch.mdate2jdate(t);
            matlabTimeString = datestr(t,irisFetch.DATE_FORMATTER);
            criteria = edu.iris.dmc.ws.criteria.RespCriteria();            
            criteria.setEndTime(javadateTime);
            reconvertedMatlabTime = ...
               datestr(irisFetch.jdate2mdate(javadateTime),irisFetch.DATE_FORMATTER);
            %urlString = char(criteria.toUrlParams().get(0));
            urlString = criteria.toUrlParams().get(0).toCharArray()';
            if ~(all(reconvertedMatlabTime == matlabTimeString));
               disp(s-fix(s));
               if datenum(reconvertedMatlabTime) > datenum(t)
                  fprintf('^ ');
               else
                  fprintf('v ');
               end
            fprintf('InputTime: %s  ; jDateTime: %s ; millis: %d\nReConvert: %s\nURL: %s\n',...
               matlabTimeString, ...
               javadateTime.toGMTString.toCharArray()', ...
               rem(javadateTime.getTime(),1000),...
               reconvertedMatlabTime,...               
               urlString);
            end
         end
      end
      
      
   end % static methods
   
   %%
   methods(Static, Access=protected)
      
      
      
      function myDateStr = makeDateStr(dateInput)
         if isnumeric(dateInput)
            myDateStr = datestr(dateInput, irisFetch.DATE_FORMATTER);
         elseif ischar(dateInput)
            myDateStr = dateInput;
         end
      end
      
      function d = jArrayList2complex(jArrayList)
         % for use on ArrayList objects containing things with getReal() and getImaginary()
         %  edu.iris.dmc.ws.sacpz.model.Pole
         %  edu.iris.dmc.ws.sacpz.model.Zero
         %  edu.iris.dmc.ws.station.model.ComplexNumber
         r= zeros(jArrayList.size(),2);
         for n=1:jArrayList.size()
            r(n,:)=double([jArrayList.get(n-1).getReal(), jArrayList.get(n-1).getImaginary()]);
         end
         
         if any(r(:,2))
            d = complex(r(:,1),r(:,2));
         else
            d = r(:,1);
         end
      end
      
      function mts = convertTraces(traces)
         %irisFetch.convertTraces converts traces from java to a matlab structure
         %USAGE:
         %  mts = convertTraces(traces) where TRACES a java trace
         %  class.
         
         blankSacPZ = struct('units','','constant',[],'poles',[],'zeros',[]);
         
         blankTrace = struct('network','','station','','location',''...
            ,'channel','','quality','',...
            'latitude',0,'longitude',0,'elevation',0,'depth',0,...
            'azimuth',0,'dip',0,...
            'sensitivity',0,'sensitivityFrequency',0,...
            'instrument','','sensitivityUnits','UNK',...
            'data',[],'sampleCount',0,'sampleRate',nan,...
            'startTime',0,'endTime',0,'sacpz',blankSacPZ);
         mts=blankTrace;
         for i = 1:length(traces)
            mt=blankTrace;
            mt.network  = char(traces(i).getNetwork());
            mt.station  = char(traces(i).getStation());
            mt.location = char(traces(i).getLocation());
            mt.channel  = char(traces(i).getChannel());
            
            mt.quality  = char(traces(i).getQuality());
            
            mt.latitude  = traces(i).getLatitude();
            mt.longitude = traces(i).getLongitude();
            mt.elevation = traces(i).getElevation();
            mt.depth     = traces(i).getDepth();
            mt.azimuth   = traces(i).getAzimuth();
            mt.dip       = traces(i).getDip();
            
            mt.sensitivity = traces(i).getSensitivity();
            mt.sensitivityFrequency = traces(i).getSensitivityFrequency();
            
            mt.instrument = traces(i).getInstrument().toCharArray()';
            mt.sensitivityUnits = traces(i).getSensitivityUnits().toCharArray()';
            mt.data = traces(i).getData();
            
            mt.sampleCount = traces(i).getSampleCount();
            mt.sampleRate = traces(i).getSampleRate();
            
            startDateString = char(traces(i).getStartTime().toString());
            endDateString = char(traces(i).getEndTime().toString());
            
            mt.startTime = datenum(startDateString, irisFetch.DATE_FORMATTER);
            mt.endTime = datenum(endDateString, irisFetch.DATE_FORMATTER);
            
            try
               jsacpz = traces(i).getSacpz();
            catch er
               if strcmp(er.identifier,'MATLAB:noSuchMethodOrField')
                  warning('IRISFETCH:convertTraces:noGetSacPZmethod',...
                     'probably using older verision of the ws-library. please retrieve the latest version');
                  jsacpz = [];
               else
                  rethrow(er)
               end
            end
            if ~isempty(jsacpz)
               sacpz.units = char(traces(i).getSacpz().getInputUnit());
               sacpz.constant = traces(i).getSacpz().getConstant();
               sacpz.poles= irisFetch.jArrayList2complex(traces(i).getSacpz().getPoles());
               sacpz.zeros= irisFetch.jArrayList2complex(traces(i).getSacpz().getZeros());
               mt.sacpz = sacpz;
            end
            mts(i) = mt;
         end
      end
      
      
      %----------------------------------------------------------------
      % DATE conversion routines
      %
      % 1970-01-01 is datenum 719529; there are 86400000 ms in a day.
      %
      % Java classes that can be used:
      %     java.sql.Timestamp : handles down to nanosecond
      %     java.util.Date     : handles milliseconds
      %
      % MATLAB itself is only accurate to the .01 milliseconds
      %----------------------------------------------------------------
      
      function javadate = mdate2jdate(matlabdate)
         %mdate2jdate converts a matlab date to a java Date class
         % TRUNCATES TO Milliseconds
         
         % 10 Oct 2012 
         % changed to use Calendar java class with milliseconds because the wrong date
         % (apparently off by 1 second) was being created
         
         if ischar(matlabdate)
            matlabdate = datenum(matlabdate);
         end
         if ~isnumeric(matlabdate) || ~isscalar(matlabdate)
            error('IRISFETCH:mdate2jdate:incorrectDateFormat',...
               'A scalar matlab datenum was expected, but a different kind of value was received.');
         end
        
         jmillis = ((matlabdate-irisFetch.BASE_DATENUM) * irisFetch.MS_IN_DAY) + .5 ; % add 0.5 to keep it in sync.
     
         %timestamp = java.sql.Timestamp(jmillis);
         javadate = java.util.Date(jmillis); %convert to a Date, loosing nanosecond precision
      end
      
      function matlabdate = jdate2mdate(javadate)
         persistent formatter
         if isempty(formatter)
            formatter = java.text.SimpleDateFormat('yyyy-MM-dd HH:mm:ss.SSS');
         end
         % jdate2mdate converts a java Date class to a matlab datenum
         % NOTE: Matlab cannot provide nanosecond resolution, though it can come close...
         % by
         if isa(javadate,'java.sql.Timestamp') %nanosecond precision
            matlabdate =  datenum([1970 1 1 0 0 (fix(javadate.getTime()/1000) + javadate.getNanos / 1000000000) ]);
            %matlab dates do not have nanosecond accuracy
         elseif isa(javadate,'java.util.Date') %millisecond precision
            matlabdate= formatter.format(javadate).toCharArray()'; % might need to transpose.
            %matlabdate= datenum([1970 1 1 0 0 (javadate.getTime()/1000)]);
         % else 
         end
%          try
%             % matlabdate = irisFetch.BASE_DATENUM + (javadate.getTime()) / irisFetch.MS_IN_DAY;
%             matlabdate=datestr(matlabdate,irisFetch.DATE_FORMATTER);
%          catch je
%             warning(je)
%             matlabdate = [];
%          end
      end
      
      %----------------------------------------------------------------
      % Look for GET / SET methods for the class.
      %----------------------------------------------------------------
      function M = getSettableFields(obj)
         % strip the first 3 letters off the field ('set')
         M = irisFetch.getSetters(obj);
         for n=1:numel(M)
            M(n) = {M{n}(4:end)};
         end
      end
      
      function [M, argType] = getSetters(obj)
         [M, argType] = irisFetch.getMethods(obj,'set');
      end
      
      function [methodList, argType] = getMethods(obj,searchPrefix)
         persistent className methodsAndArguments
         if isempty(className)
            className = {''};
            methodsAndArguments = {{''},{''}};
         end
         
         thisClass = class(obj);
         TF = strcmp(thisClass, className);
         if any(TF) % shortcut if this has been done before
            methodList = methodsAndArguments{TF,1};
            argType = methodsAndArguments{TF,2};
            return
         else
            loc = numel(className)+1;
         end
         
         
         argType = {}; %methodList = [];
         M = methods(obj);
         M2 = methods(obj,'-full');
         idx = strncmp(searchPrefix,M, length(searchPrefix));
         methodList = M(idx);
         argList = M2(idx);
         
         p1=strfind(argList,'(');
         p2=strfind(argList,')');
         for n=1:numel(argList)
            argType(n) = {argList{n}(p1{n}+1:p2{n}-1)};
         end
         
         className(loc) = {thisClass};
         methodsAndArguments(loc,1) = {methodList};
         methodsAndArguments(loc,2) = {argType};
         
      end
      %%
      %================================================================
      %----------------------------------------------------------------
      % BEGIN: PARSING ROUTINES
      %----------------------------------------------------------------
      
      function [getterList, fieldList] = getMethodsAndFields(obj)
         
         
         % this function uses a cache to speed up the retrieval of
         % get_methods and fieldnames.
         
         persistent className
         persistent methodL
         persistent fieldL
         
         if isempty(className)
            className = {''};
            methodL = {''};
            fieldL = {''};
         end
         
         thisClass = class(obj);
         
         TF = strcmp(thisClass, className);
         
         if any(TF) % shortcut if this has been done before
            getterList = methodL{TF};
            fieldList = fieldL{TF};
            return
         else
            loc = numel(className)+1;
         end
         
         allMethods = methods(obj);
         getterList = allMethods(strncmp('get',allMethods, 3));
         
         % filter classes need to have class names for them to make sense
         % to the user.
         if isa(obj,'edu.iris.dmc.ws.station.model.Filter')
            getterList = getterList(...
               ~( strcmp('get',getterList) | ...
               strcmp('getAny',getterList) ));
            
         else
            getterList = getterList(...
               ~( strcmp('get',getterList) | ...
               strcmp('getClass',getterList) | ...
               strcmp('getAny',getterList) ));
         end
         % eliminate recursions
         switch thisClass
            case 'edu.iris.dmc.ws.station.model.Station'
               n = strcmp(getterList,'getNetwork');
            case 'edu.iris.dmc.ws.station.model.StationEpoch'
               n = strcmp(getterList,'getStation');
            case 'edu.iris.dmc.ws.station.model.Channel'
               n = strcmp(getterList,'getStationEpoch');
            case 'edu.iris.dmc.ws.station.model.ChannelEpoch'
               n = strcmp(getterList,'getChannel');
            case 'edu.iris.dmc.ws.station.model.Response'
               n = strcmp(getterList,'getChannelEpoch');
            case 'edu.iris.dmc.ws.station.model.Sensor'
               n = strcmp(getterList,'getChannelEpoch');
            otherwise
               n=[];
         end
         getterList(n) = [];
         
         fieldList = strrep(getterList,'get',''); %get rid of 'get'
         
         className(loc) = {thisClass};
         methodL(loc) = {getterList};
         fieldL(loc) = {fieldList};
      end
      
      function myStruct = parseObjectViaGetMethods(thisObj)
         % This routine should only be called for objects. Not for arrays
         % NOTE: assumes a single/scalar object.
         myStruct=[];
         [getterList, fieldnameList] = irisFetch.getMethodsAndFields(thisObj);
         
         for idx = 1 : numel(getterList)
            value = thisObj.(getterList{idx});
            if isjava(value)
               myStruct.(fieldnameList{idx}) = irisFetch.parse(value);
            else
               myStruct.(fieldnameList{idx}) = value;
            end
         end
         %
         % disp(myStruct);
         %
         
      end %function_parseObjectViaGetMethods
      
      
      function myGuts = parse(obj)
         % parse this object, first dealing with java built-in classes,
         % then proceeding to the iris classes
         
         myClass = class(obj);
         % disp(myClass); %DEBUG CODE!!!!!!
         firstbit = myClass(1:find(myClass == '.',1,'first')-1);
         %switch strtok(myClass,'.')
         switch firstbit
            case 'java' % deal with a built-in java class
               switch myClass
                  case {'java.lang.String'}
                     %myGuts = char(obj);
                     myGuts = obj.toCharArray()';
                  case 'java.lang.Class'
                     myGuts = obj.getCanonicalName.toCharArray()';
                     
                  case {'java.lang.Double',...
                        'java.lang.Integer',...
                        'java.math.BigInteger',...
                        'java.lang.Long'}
                     myGuts = obj.doubleValue;
                     
                  case 'java.util.ArrayList' % this was an array of arrays
                     % myGuts = irisFetch.parseCollection(obj);
                     % myGuts = irisFetch.parseArrayList(obj, level);
                     
                     if obj.isEmpty();
                        myGuts = [];
                     else
                        
                        
                        switch class(obj.get(0))
                           
                           case {'edu.iris.dmc.ws.sacpz.model.Pole', ...
                                 'edu.iris.dmc.ws.sacpz.model.Zero', ...
                                 'edu.iris.dmc.ws.station.model.ComplexNumber'}
                              myGuts = irisFetch.jArrayList2complex(obj);
                              
                           case {'java.lang.Double',...
                                 'java.lang.Integer',...
                                 'java.math.BigInteger',...
                                 'java.lang.Long',...
                                 'double'}
                              myGuts = double(obj.toArray(javaArray('java.lang.Double',obj.size())));
                              
                           otherwise
                              for n = obj.size:-1:1
                                 mG = irisFetch.parse(obj.get(n-1));
                                 try
                                    myGuts(n) = mG;
                                 catch er
                                    % differing versions of matlab have different spellings of this error.
                                    switch er.identifier
                                       case {'MATLAB:heterogeneousStrucAssignment', 'MATLAB:heterogenousStrucAssignment'}
                                          f = fieldnames(mG);
                                          for z=1:numel(f)
                                             myGuts(n).(f{z}) = mG.(f{z});
                                          end
                                       otherwise
                                          rethrow (er)
                                    end %switch
                                 end %try/catch
                                 
                              end %obj loop
                        end
                     end %endif obj is empty
                     
                  case 'java.util.Date'
                     myGuts = irisFetch.jdate2mdate(obj);
                     
                  otherwise
                     disp(['not sure how to deal with JAVA class :', myClass]);
               end
               
            case 'edu' % deal with one of our classes
               switch myClass
                  case {'edu.iris.dmc.ws.station.model.Sensitivity','edu.iris.dmc.ws.station.model.Sensor'}
                     
                     % in versions of irisFetch prior to 1.2, these were
                     % add these particular classes to the existing structure
                     % WITHOUT nesting to another level. Instead, use the
                     % class name as the prefix.
                     
                     % beware of recursions!!!
                     [getterList, fieldnameList] = irisFetch.getMethodsAndFields(obj);
                     prefix = myClass( find( myClass=='.', 1, 'last') + 1:end);
                     fieldnameList = strcat(prefix, fieldnameList);
                     
                     for idx = 1 : numel(getterList)
                        myGuts.(fieldnameList{idx}) = irisFetch.parse(obj.(getterList{idx}));
                     end
                  case 'edu.iris.dmc.ws.event.model.Arrival'
                     % in the hopes of speeding this up, tackle it directly
                     % before adding this, a few days in feb 2010 minmag 5 takes a long
                     % time.
                     myGuts.Distance     = obj.getDistance.doubleValue();
                     myGuts.TimeResidual = obj.getTimeResidual.doubleValue(); 
                     % myGuts.Phase        = char(obj.getPhase);
                     myGuts.Phase2        = obj.getPhase.toCharArray()';
                     myGuts.Azimuth      = obj.getAzimuth.doubleValue();
                     myGuts.Picks        = obj.getPicks;
                     if myGuts.Picks.isEmpty()
                        myGuts.Picks=[];
                     else
                        error('Not ready for this one yet')
                     end
                     myGuts.PublicId     = obj.getPublicId.toCharArray()';
                     
                  otherwise
                     myGuts = irisFetch.parseObjectViaGetMethods(obj);
               end
            otherwise
               % take best guess
               myGuts = irisFetch.parseObjectViaGetMethods(obj);
         end
         
      end % fuction_parse
      
      function myguts = parseArrayList(obj)
         % it is known that obj is of class java.util.ArrayList
         if obj.isEmpty(); myguts = []; return; end;
         switch class(obj.get(0))
            
            case {'edu.iris.dmc.ws.sacpz.model.Pole', ...
                  'edu.iris.dmc.ws.sacpz.model.Zero', ...
                  'edu.iris.dmc.ws.station.model.ComplexNumber'}
               myguts = irisFetch.jArrayList2complex(obj);
               
            case {'java.lang.Double',...
                  'java.lang.Integer',...
                  'java.math.BigInteger',...
                  'java.lang.Long',...
                  'double'}
               myguts = double(obj.toArray(javaArray('java.lang.Double',obj.size())));

            otherwise
               for n = obj.size:-1:1
                  myguts(n) = irisFetch.parse(obj.get(n-1));
               end
         end
         
      end
      
      %----------------------------------------------------------------
      % END: PARSING ROUTINES
      %----------------------------------------------------------------
      %================================================================
      %%
      
      function criteria = addCriteria(criteria, value, addMethod)
         %used to add Sta, Net, Loc, Chan to criteria
         % For example:
         %   singleAddToCriteria(criteria, {'OKCF','MNO?'},'addStation')
         % will invoke:
         %   criteria.addStation('OKCF').addStation('MNO?')
         if isempty(value)
            return %do nothing
         end
         if ~iscell(value)
            % it's probably a string
            criteria.(addMethod)(value);
         else
            % a cell may have multiple values, use 'em all.
            for n=1:numel(value)
               criteria.(addMethod)(value{n})
            end
         end
      end
      
      
      function criteria = setBoxCoordinates(criteria, thisValue)
         % setBoxCoordinates (minLat, maxLat, minLon, maxLon)
         % values of 'NAN' are ignored
         if numel(thisValue) ~=4
            error('IRISFETCH:setBoxCoordinates:InvalidParameterCount',...
               'Expected [minLat, maxLat, minLon, maxLon]');
         end
         setMethods = {'setMinimumLatitude','setMaximumLatitude',...
            'setMinimumLongitude','setMaximumLongitude'};
         for n=1:numel(setMethods)
            if ~isnan(thisValue(n))
               criteria.(setMethods{n})(java.lang.Double(thisValue(n)));
            end
         end
      end
      
      function criteria = setCriteria(criteria, paramList)
         
         %----------------------------------------------------------
         % The following code allows for open-ended search criteria
         % allowing it to change whenever the java library changes
         %
         % Instead of hard-coding each Setter, I query the java class to
         % find out its methods, then keep the ones that start with
         % "set".
         %
         % I also find out what the input parameters are for each, and
         % use that to properly create/cast the data.  Without doing
         % this, neither dates nor doubles would work. Boo.
         %----------------------------------------------------------
         
         % Get a list of parameters, their set functions, and input
         % types, and do it outside the loop so they are not needlessly
         % rerun
         
         [allSetMethods, argType] = irisFetch.getSetters(criteria);
         settableFieldnames = irisFetch.getSettableFields(criteria);
         allMethods = methods(criteria);
         
         while ~isempty(paramList) && numel(paramList) >= 2
            
            % Grab the parameter pair, then remove from parameter list
            thisParam = paramList{1};
            thisValue = paramList{2};
            paramList(1:2)=[];
            
            indexOfMethod = strcmpi(thisParam,settableFieldnames);
            if any(indexOfMethod)
               
               setMethod = allSetMethods{indexOfMethod};
               switch argType{indexOfMethod}
                  
                  case 'java.util.Date'
                     criteria.(setMethod)(irisFetch.mdate2jdate(thisValue));
                     criteria.toUrlParams;
                  case 'java.lang.Double'
                     criteria.(setMethod)(java.lang.Double(thisValue));
                     
                  case 'java.lang.Long'
                     criteria.(setMethod)(java.lang.Long(thisValue));
                     
                  case 'java.lang.Integer'
                     criteria.(setMethod)(java.lang.Integer(thisValue));
                     
                  otherwise
                     disp('Unanticipated argument type... trying');
                     criteria.(setMethod)(thisValue);
               end
               continue;
            end
            % we shall only pass this point if existing methods were not used
            
            switch lower(thisParam)
               %handle special cases
               case 'boxcoordinates'
                  criteria = irisFetch.setBoxCoordinates(criteria, thisValue);
                  
               case 'radialcoordinates'
                  criteria.setLatitude(java.lang.Double(thisValue(1)));
                  criteria.setLongitude(java.lang.Double(thisValue(2)));
                  criteria.setMaximumRadius(java.lang.Double(thisValue(3)));
                  if numel(thisValue) ==4 && ~isnan(thisValue(4))
                     criteria.setMinimumRadius(java.lang.Double(thisValue(4)));
                  end
                  
               otherwise
                  % this will blow up if java doesn't recongize
                  % thisValue
                  criteria.(allMethods{strcmpi(allMethods,thisParam)})(thisValue);
            end
         end
      end
      
      %--------------------------------------------------------------------
      
   end %static protected methods
end
