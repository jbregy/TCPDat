% TCPDat: Tropical Cyclone Precipitation Dataset
% Public version 1.0 release 1.0
% Date updated: 03/02/2020 (dd/mm/yyyy)

% If you use this code or the pre-generated TCPDat data, please include the
% following citation (or whichever citation style your journal of choice requies):

% Bregy, J.C., J.T. Maxwell, S.M. Robeson, J.T. Ortegren, P.T. Soulé, and
% P.A. Knapp, 2019: Spatiotemporal variability of tropical cyclone
% precipitation using a high-resolution, gridded (0.25° x 0.25°) dataset
% for the eastern United States, 1948–2015. Journal of Climate, 
% 33(5), 1803–1819, doi:10.1175/JCLI-D-18-0885.1.

% Stable link to the original paper: https://doi.org/10.1175/JCLI-D-18-0885.1


% This code generates a gridded tropical cyclone precipitation (TCP)
% product for the eastern United States (i.e., east of the Rocky Mountain
% continental divide) from 1948 to present-1 (or present-2, depending on
% the DOY). This script requires two (2) publically available datasets in
% order to create TCPDat. The first is HURDAT2, specifically the North
% Atlantic version of IBTrACS (International Best Track Archive for Climate
% Stewardship). The script is designed to work with the netCDF version of
% IBTrACS. You can find the latest version of IBTrACS at:
% Cite IBTrACS accordingly.
% The other dataset required is the Climate Prediction Center Unified
% Guage-Based Analysis of Daily Precipitation over CONUS (CPC US Unified)
% from the Physical Sciences Division of the Earth Systems Research
% Laboratory. The latest version can be found here: 
% https://www.esrl.noaa.gov/psd/thredds/catalog/Datasets/cpc_us_precip/catalog.html
% (Note: the years 2007 onward are located in the RT subfile. CPC started doing something different.
% It's inconvenient, but such is life.) Like IBTrACS, CPC US Unified should be a netCDF. You will have to
% download gridded precipitation for each year that you are interested in.
% The script is written to cover every year starting at 1948 (per CPC US
% Unified). In future iterations, I plan to have query statements that
% allow you to define your time range, but plans always look great on paper.
% Anyway, cite CPC US Unified accordingly.

% This is my attempt at annotating this code. It is going to be rough; this
% is literally the first time I have ever done something like this. Let the
% GitHub gods have mercy upon my digital soul. If there are any questions
% about anything that is in the code, or if you have any suggestions to
% improve the code, please do not hesitate to get in touch with me. My most
% frequently checked email is jbregy@indiana.edu, but I also have a stable
% address at joshua.bregy@gmail.com. I suppose you could get in touch via
% GitHub, ResearchGate, or Twitter (@prehistormic) if you would like.

%% Read this before running the script

% You will need to create several folders for this script to work. Each
% folder houses different data, and the code moves between different folders 
% to either use or save datasets. I wrote it this way to force myself to be
% organized (normally, I'm not, but I knew my future self would genuinely
% appreciate it). Here are the folders that you will need:

% 1) A folder containing your CPC US Unified files. I called my
% CONUS_daily_precip_1948_2016. (I accidentally downloaded p2016.nc, but
% HURDAT2 had not been updated to include 2016.) I will be changing the
% file name in the command to CONUS_daily_precip, so I suggest that you
% create a file with the same name. That, or change the name in the code
% yourself. Alternatively, you could just comment out directory changes
% that are used to access data as long as the folder containing your
% precipitation (and I suppose IBTrACS, too) is included whenever you set
% your path.

% 2) A folder for the TCP data. I have a line in the code that creates that
% folder. It is named TCP_data. 
%% TC data extraction

clear % Do what you want with this. It's a force of habit by now.
close % Same thing.

first_run = input('Do you already have a folder named TCP_Data in this parent folder or the in some path? 1 = Yes; 0 = No (sorry, I am lazy)')
if first_run == 0
    mkdir TCP_Data
else
end
%ncdisp('Basin.NA.ibtracs_wmo.v03r09.nc') % I usually just have this in so
% I can see what variables I'm working with. 
TCfile_name = input('What is the name of your netCDF file?','s') %Come back 
% to this in order to incoporate and streamline it.
filename_ext = [TCfile_name,'.nc'];
ncid = netcdf.open(filename_ext,'NC_NOWRITE'); % Simply calling the file
% name. I don't want to type it repeatedly. Change as you see fit. Also,
% this should open the user specified netCDF file that the user is working 
% with. If not, look at the example I have on Line 60 (or the line below). 
% Just change the file name if you are going to use that line of code 
% instead. Also, comment out the ncid on Line 53.

%ncid = netcdf.open('Basin.NA.ibtracs_wmo.v03r09.nc','NC_NOWRITE'); % Same
% thing, but done in a less eloquent manner. 
varid = netcdf.inqVarID(ncid,'season');
season = double(netcdf.getVar(ncid,varid)); % Year based on season (01 January
% for the Northern Hemisphere). Yes, TCs happen outside of JJASON.
% season = double(netcdf.getVar(ncid,3)); % Remnant call.
st_yr = input('What year do you want to start at?') % Just adding this in 
% the event that you want to start at a different year. Remember, TCPDat
% can only go back to 1948 because that is the oldest year in CPC US
% Unified.
start_season = find(season == st_yr,1,'first'); % Index location of 1948. If
% this doesn't work, try the commented out line below.
%s48 = find(season == 1948,1, 'first'); % Index location of 1948.
ed_yr = input('What year do you want to end at?') % Same thing as st_yr.
end_season = find(season == ed_yr,1,'last'); % Index location of the last year.
% If this doesn't work, then try the commented out line below.
%e15 = find(season == 2015,1,'last'); % Index location of 2015. Change as
% you see fit.
year = season(start_season:end_season,:); % Extracts the years of interested. For
% example, if you wanted to look at TCP from 1948 to 2015, then that's what
% would be extracted. If that does not work, try the commented out line
% below.
%year = season(s48:e15,:); % 1948 to 2015. Use if the other approach
%doesn't work.

%varid = netcdf.inqVarID(ncid,'nature'); 
%nature = double(netcdf.getVar(ncid,varid)); % Storm nature: 0 = TS - tropical (i.e. hurricane,
%tropical cyclone,tropical storm, tropical depression, typhoon, and super
%typhoon) -127 = NaN
% nature = double(netcdf.getVar(ncid,17)); % Remnant call.
%nature = nature(:,start_season:end_season); % Gets the TCs.
%nature= nature(:,s48:e15);
%nature(nature == -127) = NaN; % Some storms have shorter lifespans than others, or
% the instruments might have failed. Either way, missing data (-127) is now
% listed as NaN.

% storm_sn = netcdf.getVar(ncid,0); % Storm serial number. I did not use
% this, but if you need it, you have it. The serial number is made of
% characters.
% storm_sn = storm_sn(:,start_season:end_season); % Same thing. If this
% doesn't work, then use the one on the line below, provided you've
% switched to s48 and e15 above.
% storm_sn = storm_sn(:,s48:e15); % Same thing.

%varid = netcdf.inqVarID(ncid,'name');
%name = netcdf.getVar(ncid,varid); % Storm name. All storms were unnamed prior to
% 1950. Nice to have, but not necessary (date confirmation). Characaters.
% name = netcdf.getVar(ncid,1); % Remnant call.
%name = name(:,start_season:end_season); % Just isolating the storm names. 
% Again, not necessary, but can be useful for date confirmation. Switch to
% the commented out line below if it doesn't work.
%name = name(:,s48:e15);

%varid = netcdf.inqVarID(ncid,'numobs'); 
%numObs = double(netcdf.getVar(ncid,varid)); % Number of observations for the 
% storm. I suppose it is useful for the spline.
% numObs = double(netcdf.getVar(ncid,2)); % Remnant call.
%numObs = numObs(start_season:end_season,:); % You know the drill. If it
% doesn't work, then just use the commented out one below.
% numObs = numObs(s48:e15,:);
%varid = netcdf.inqVarID(ncid,'track_type');
%track_type = double(netcdf.getVar(ncid,varid)); % 0 = main (cyclogenesis (cg)
% to cyclolysis (cl)), 1 = merge (cg to merge), 2 = split (split to cl),
% 3 = other (split to merger). Might not be necessary.
% track_type = double(netcdf.getVar(ncid,4)); % Remnant call.
%track_type = track_type(start_season:end_season,:); % Same thing. Just switch
% it if it does not work. 
% track_type = track_type(s48:e15,:);

%genesis_basin = double(netcdf.getVar(ncid,5)); % 0 = North Atlantic (NA), 
% 12 = Caribbean Sea (CS),13 = Gulf of Mexico (GM). Useful if you want to
% do some mapping or whatever.
%genesis_basin = genesis_basin(start_season:end_season,:); % Switch if it
% doesn't work.
% genesis_basin = genesis_basin(s48:e15,:);

%num_basins = double(netcdf.getVar(ncid,6)); % Number of basins through which
% the storm passes. Not really necessary.
%num_basins = num_basins(start_season:end_season,:); % Same thing. Switch.
% num_basins = num_basins(s48:e15,:);
%varid = netcdf.inqVarID(ncid,'basin');
%basin = double(netcdf.getVar(ncid,varid)); % Based on present location. 0 = NA,
% 12 = CS, 13 = GM -127 = NaN
% basin = double(netcdf.getVar(ncid,7)); % Remnant call.
%basin = basin(:,start_season:end_season); % Getting the specific basins.
% Switch to below if necessary.
%basin = basin(:,s48:e15);
%basin(basin == -127) = NaN; % Missing values (-127) converted to NaNs.

%wind_avg_period = netcdf.getVar(ncid,8); % Wind speed averaging period. 
% Units: min. -127 = NaN

%source = netcdf.getVar(ncid,9); %Source name. Matches dimensions in source_* variables
varid = netcdf.inqVarID(ncid,'time');
time_wmo = double(netcdf.getVar(ncid,varid)); % Modified Julian Day. 
% Units: days since 1858-11-17 00:00:00. ~9,97*10^36 = NaN
% time_wmo = double(netcdf.getVar(ncid,10)); % Remnant call.
time_wmo = time_wmo(:,start_season:end_season);
%time_wmo = time_wmo(:,s48:e15);
time_wmo(time_wmo < -9900) = NaN;%== 9.969209999999999e+36) = NaN;

varid = netcdf.inqVarID(ncid,'lat');
lat_wmo = double(netcdf.getVar(ncid,varid)); % Storm center latitude. 
% Units: deg_N. scale factor: 0.01 -32767 = NaN
% lat_wmo = double(netcdf.getVar(ncid,11)); % Remnant call.
lat_wmo = lat_wmo(:,start_season:end_season);
%lat_wmo = lat_wmo(:,s48:e15);
lat_wmo(lat_wmo == -9999) = NaN;%-32767) = NaN;
%lat_wmo = lat_wmo*0.01; % Scaling (0.1) per the notes in IBTrACS.
% I don't think I need the scaling factor

varid = netcdf.inqVarID(ncid,'lon');
lon_wmo = double(netcdf.getVar(ncid,varid)); % Storm center longitude. 
% Units: deg_E. scale factor: 0.01 -32767 = NaN
% lon_wmo = double(netcdf.getVar(ncid,12)); % Remnant call.
lon_wmo = lon_wmo(:,start_season:end_season);
%lon_wmo = lon_wmo(:,s48:e15);
lon_wmo(lon_wmo == -9999) = NaN;%-32767) = NaN;
%lon_wmo = lon_wmo*0.01; % Same thing as latitude.
%I don't think I need the scaling factor.

varid = netcdf.inqVarID(ncid,'wmo_wind');
wind_wmo = double(netcdf.getVar(ncid,varid)); % max sustained wind 
% Units: kt scaling factor: 0.1 -32767 = NaN
% wind_wmo = double(netcdf.getVar(ncid,14)); % Remnant call.
wind_wmo = wind_wmo(:,start_season:end_season);
%wind_wmo = wind_wmo(:,s48:e15);
wind_wmo(wind_wmo == -9999) = NaN;%-32767) = NaN;
wind_wmo = wind_wmo*0.51; % converting to m/s
%wind_wmo = (wind_wmo*0.1)*0.51; % scaling (*0.1) and converting to m/s (*0.51).
% Wind scaling might not be needed. It's not clear in the updated IBTrACS.

varid = netcdf.inqVarID(ncid,'wmo_pres');
pres_wmo = double(netcdf.getVar(ncid,varid)); % Minimum central pressure. 
% Units: mb scaling factor: 0.1 -32767 = NaN
% pres_wmo = double(netcdf.getVar(ncid,15)); % Remnant call.
pres_wmo = pres_wmo(:,start_season:end_season);
%pres_wmo = pres_wmo(:,s48:e15);
pres_wmo(pres_wmo == -9999) = NaN;%-32767) = NaN;
%pres_wmo = pres_wmo*0.1; % Scaling (0.1)
% Scaling might not be needed. It's not clear in the updated IBTrACS.

%varid = netcdf.inqVarID(ncid,'subbasin');
%sub_basin = double(netcdf.getVar(ncid,varid)); % Same classification as the 
% other basin variabiles. -127 = NaN
% sub_basin = double(netcdf.getVar(ncid,16)); % Remnant call.
%sub_basin = sub_basin(:,start_season:end_season);
%sub_basin = sub_basin(:,s48:e15);
%sub_basin(sub_basin == -127) = NaN;

%source_wmo = double(netcdf.getVar(ncid,18)); % data source (updated HURDAT
% from NHC and HRD) -127 = NaN. Honestly, I'm not sure what this is really
% even supposed to be.
%source_wmo = source_wmo(:,start_season:end_season);
%source_wmo = source_wmo(:,s48:e15);
%source_wmo(source_wmo == -127) = NaN;

varid = netcdf.inqVarID(ncid,'dist2land');
dist2land = double(netcdf.getVar(ncid,varid)); % Distance to land 
% Units: km -999 = NaN
% dist2land = double(netcdf.getVar(ncid,19)); % Remnant call.
dist2land = dist2land(:,start_season:end_season);
%dist2land = dist2land(:,s48:e15);
dist2land(dist2land == -9999) = NaN;

varid = netcdf.inqVarID(ncid,'landfall');
landfall = double(netcdf.getVar(ncid,varid)); % Minimum distance to land until
% next report (0=landfall) units: km -999 = NaN
% landfall = double(netcdf.getVar(ncid,20)); % Remnant call.
landfall = landfall(:,start_season:end_season);
%landfall = landfall(:,s48:e15);
landfall(landfall == -9999) = NaN;

varid = netcdf.inqVarID(ncid,'usa_sshs');
nature_num = double(netcdf.getVar(ncid,varid)); % Using this instead of nature because it's a numerical representation.
nature_num = nature_num(:,start_season:end_season);
nature_num(nature_num == -15) = NaN;

%% Isolate tropical systems

nontc_loc = find(nature_num == -5 | nature_num == -4 | nature_num == -3 | nature_num == -2); % Use this instead of remnant call below.
%nontc_loc = find(nature == 1 | nature == 2 | nature == 3 | nature == 4 |...
 %   nature == 5 | nature == 6 | nature == 7); % identify the location of nonTCs
% nature(nontc_loc) = NaN; you can apply this to any of the variables,
% which would eliminate the need for the for loop below. For now, let's
% keep it until you're making the small changes.
[natx,naty] = size(nature_num);%nature);
natsize = [natx,naty];
[nontc_row_loc,nontc_col_loc] = ind2sub(natsize,nontc_loc); %convert linear 
% indices to subscript indices

for i = 1:naty % replace with column variable for nature
    for j = 1:natx % same thing, but with rows
        if (nature_num(j,i) < -1) && (nature_num(j,i) > -6) % This just focuses on the 
            % nonTC entries. Recall that 0 represents tropical systems.
            % Once nonTC entries are found (including extratropical or
            % subtropical stages of a TC), they are populated with NaNs 
            % (see below). The reason why this is in a loop is because we're 
            % going through each storm and its associated entries.
        % if (nature(j,i) > 0) && (nature(j,i) < 8) % Remnant call.
            %nature(j,i) = NaN;
            nature_num(j,i) = NaN;
            %basin(j,i) = NaN;
            dist2land(j,i) = NaN;
            landfall(j,i) = NaN;
            lat_wmo(j,i) = NaN;
            lon_wmo(j,i) = NaN;
            pres_wmo(j,i) = NaN;
            %sub_basin(j,i) = NaN;
            time_wmo(j,i) = NaN;
            wind_wmo(j,i) = NaN;
        end
    end
end

time_wmo_conv = datetime(time_wmo,'ConvertFrom','modifiedjuliandate'); % MJD 
% to Gregorian and UTC. Time variables are extracted below.
year_utc = time_wmo_conv.Year; month_utc = time_wmo_conv.Month;
day_utc = time_wmo_conv.Day; hour_utc = (time_wmo_conv.Hour)*100;


% ID storms that are 223 km or less from any land. This distance is based
% on the average radius of a TC rainfield per Matyas (2010). See Bregy et al.,
% 2019 for discussion of this rainfield. Although we used 223 km to build the 
% dataset, you are free to change the radius as you see fit. 
count = 0; % this and the for-if loop determines the number of storms within 
% the specified distance from land.
user_distance = input('What rainfield radius would you like to use?')
for i = 1:naty
    if any(dist2land(:,i) <= user_distance) % This should be able to use your 
        % input but life likes to be inconvient. If it doesn't work, use the 
        % commented out line of code in the if statement (immediately below).
    %if any(dist2land(:,i) <= 223)
        count  = count + 1;
    else
        count = count + 0;
    end
end

% The variables (tc*) below were created to build TCPDat.
tcdist2land = zeros(j,count);
tcnature = zeros(j,count);
%tcbasin = zeros(j,count);
tclandfall = zeros(j,count);
tclat = zeros(j,count);
tclon = zeros(j,count);
tcpress = zeros(j,count);
%tcsub_basin = zeros(j,count);
tcyear = zeros(j,count);
tcmonth = zeros(j,count);
tcday = zeros(j,count);
tchour = zeros(j,count);
tcwind = zeros(j,count);
count = 1; % Used in the for loop/if statement.
for i = 1:naty
    if any(dist2land(:,i) <= user_distance) % Basically, this is identifying whether 
        % any entries for a storm is with 223 km of land. This does not 
        % distinguish between CONUS and other countries. However, CPC US Unified 
        % only covers CONUS, using NaNs elsewhere (other countries and over the ocean),
        % and therefore, we are able to extract precipitation only over CONUS.
        % Use old if statement condition if new one doesn't work: any(dist2land(:,i) <= 223)
        tcdist2land(:,count) = dist2land(:,i);
        tcnature(:,count) = nature_num(:,i);%nature(:,i);
       % tcbasin(:,count) = basin(:,i);
        tclandfall(:,count) = landfall(:,i);
        tclat(:,count) = lat_wmo(:,i);
        tclon(:,count) = lon_wmo(:,i);
        tcpress(:,count) = pres_wmo(:,i);
        %tcsub_basin(:,count) = sub_basin(:,i);
        tcyear(:,count) = year_utc(:,i);
        tcmonth(:,count) = month_utc(:,i);
        tcday(:,count) = day_utc(:,i);
        tchour(:,count) = hour_utc(:,i);
        count = count+1;
    end
end


yr_uniq = unique(tcyear(1,:));
yr_uniq = yr_uniq(1,find(yr_uniq == st_yr,1,'first'):...
    find(yr_uniq == ed_yr,1,'first')); % If this doesn't work, then use 
    % the commented out part of the code below. Change the years accordingly.
%yr_uniq = yr_uniq(1,find(yr_uniq == 1948,1,'first'):find(yr_uniq == 2015,1,'first'));
 

%% Creating TCPDat
cd CONUS_daily_precip_1948_2016\ % This line is user specific. 
% cd CONUS_daily_precip\ 
pd = dir('p*.nc');
pdl = length(pd);
cd ../ 
TCP_annual = zeros(300,120,length(yr_uniq)); % This is creating an annual sum of TCP. 
% The third dimension is based on the user input for years.
for b = 1:pdl% 58 is 2005 youcould also use the length of yr_uniq.
    
    % The if statement below is indexing everything.
    if yr_uniq(b) == 1948 % 1948 is strange and has to be singled out.
        tcyroi = tcyear(:,(find(tcyear(1,:) == yr_uniq(1,b),1,'first'):find(tcyear(1,:) == yr_uniq(1,b),1,'last'))); % Indexing TC year of interest.
        tcyroi_loc_st = find(tcyear(1,:) == yr_uniq(1,b),1,'first');
        tcyroi_loc_en = find(tcyear(1,:) == yr_uniq(1,b),1,'last');
    elseif yr_uniq(b) == 2015 | yr_uniq(b) == 2018 % 2015 is strange and has to be singled out.
        tcyroi = tcyear(:,(find(tcyear(1,:) == yr_uniq(1,b),1,'first')-1):(find(tcyear(1,:) == yr_uniq(1,b),1,'last')));
        tcyroi_loc_st = find(tcyear(1,:) == yr_uniq(1,b),1,'first')-1;
        tcyroi_loc_en = find(tcyear(1,:) == yr_uniq(1,b),1,'last');
    elseif yr_uniq(b) == 1954 % Hurricane Alice is why we can't have nice things.
        % Hurricane Alice formed on 30/12/1954 and dissipated on 6/1/1955.
        tcyroi = tcyear(:,(find(tcyear(1,:) == yr_uniq(1,b),1,'first'):(find(tcyear(1,:) == yr_uniq(1,b),1,'last')-1)));
        tcyroi_loc_st = find(tcyear(1,:) == yr_uniq(1,b),1,'first');
        tcyroi_loc_en = find(tcyear(1,:) == yr_uniq(1,b),1,'last')-1;
    else % Everything else should be good to go. There might be years that are 
        % similar to 2015, but I will update that in future version.
        tcyroi = tcyear(:,(find(tcyear(1,:) == yr_uniq(1,b),1,'first'):find(tcyear(1,:) == yr_uniq(1,b),1,'last')));
        tcyroi_loc_st = find(tcyear(1,:) == yr_uniq(1,b),1,'first');
        tcyroi_loc_en = find(tcyear(1,:) == yr_uniq(1,b),1,'last');
        nancheck_st = isnan(tcyear(1,tcyroi_loc_st-1));
        nancheck_ed = isnan(tcyear(1,tcyroi_loc_en+1));
        if nancheck_st == 1 && isequal(tcyear(1,tcyroi_loc_st),tcyear(find(tcyear(:,tcyroi_loc_st-1) == yr_uniq(1,b),1,'first'),tcyroi_loc_st-1))...
                && (nancheck_ed ~= 1 || (nancheck_ed == 1 && isequal(tcyear(1,tcyroi_loc_en),tcyear(find(tcyear(:,tcyroi_loc_en+1) == yr_uniq(1,b),1,'first'),tcyroi_loc_en+1))) == 0)
            tcyroi = tcyear(:,(find(tcyear(1,:) == yr_uniq(1,b-1),1,'last')+1):(find(tcyear(1,:) == yr_uniq(1,b),1,'last')));
            tcyroi_loc_st = find(tcyear(1,:) == yr_uniq(1,b-1),1,'last')+1;
            tcyroi_loc_en = find(tcyear(1,:) == yr_uniq(1,b),1,'last');
        elseif nancheck_ed == 1 && isequal(tcyear(1,tcyroi_loc_en),tcyear(find(tcyear(:,tcyroi_loc_en+1) == yr_uniq(1,b),1,'first'),tcyroi_loc_en+1))...
                && (nancheck_st ~= 1 || (nancheck_st == 1 && isequal(tcyear(1,tcyroi_loc_st),tcyear(find(tcyear(:,tcyroi_loc_st-1) == yr_uniq(1,b),1,'first'),tcyroi_loc_st-1))) == 0)
            tcyroi = tcyear(:,(find(tcyear(1,:) == yr_uniq(1,b),1,'first')):(find(tcyear(1,:) == yr_uniq(1,b+1),1,'first')-1));
            tcyroi_loc_st = find(tcyear(1,:) == yr_uniq(1,b),1,'first');
            tcyroi_loc_en = find(tcyear(1,:) == yr_uniq(1,b+1),1,'first')-1;
        elseif nancheck_st == 1 && nancheck_ed == 1 && isequal(tcyear(1,tcyroi_loc_st),tcyear(find(tcyear(:,tcyroi_loc_st-1) == yr_uniq(1,b),1,'first'),tcyroi_loc_st-1))...%tcyear(1,tcyroi_loc_st) == tcyear(find(tcyear(:,tcyroi_loc_st-1) == yr_uniq(1,b),1,'first'),tcyroi_loc_st-1)...
                && isequal(tcyear(1,tcyroi_loc_en),tcyear(find(tcyear(:,tcyroi_loc_en+1) == yr_uniq(1,b),1,'first'),tcyroi_loc_en+1))%tcyear(1,tcyroi_loc_en) == tcyear(find(tcyear(:,tcyroi_loc_en+1) == yr_uniq(1,b),1,'first'),tcyroi_loc_en+1);
            tcyroi = tcyear(:,(find(tcyear(1,:) == yr_uniq(1,b-1),1,'last')+1):(find(tcyear(1,:) == yr_uniq(1,b+1),1,'first')-1));
            tcyroi_loc_st = find(tcyear(1,:) == yr_uniq(1,b-1),1,'last')+1;
            tcyroi_loc_en = find(tcyear(1,:) == yr_uniq(1,b+1),1,'first')-1;
        else
        end
    end
    
    % This is extracting the data for the TCs of interest
    tcmon_oi = tcmonth(:,tcyroi_loc_st:tcyroi_loc_en);
   % tcbasin_oi = tcbasin(:,tcyroi_loc_st:tcyroi_loc_en);
    tcday_oi = tcday(:,tcyroi_loc_st:tcyroi_loc_en);
    tcdist2land_oi = tcdist2land(:,tcyroi_loc_st:tcyroi_loc_en);
    tchour_oi = tchour(:,tcyroi_loc_st:tcyroi_loc_en);
    tclandfall_oi = tclandfall(:,tcyroi_loc_st:tcyroi_loc_en);
    tclat_oi = tclat(:,tcyroi_loc_st:tcyroi_loc_en);
    tclon_oi = tclon(:,tcyroi_loc_st:tcyroi_loc_en);
    tcnature_oi = tcnature(:,tcyroi_loc_st:tcyroi_loc_en);
    tcdate_doy = day(datetime(tcyroi,tcmon_oi,tcday_oi),'dayofyear');
  
    
    % Opening URD files and extracting data.
    cd CONUS_daily_precip_1948_2016\ 
    pncid = netcdf.open(pd(b).name,'NC_NOWRITE');
    
    plat = double(netcdf.getVar(pncid,0)); %latitude of precipitation grid cells
    plon = double(netcdf.getVar(pncid,1)); %longitude of precipitation grid cells
    plon = plon-360;
    ptime = double(netcdf.getVar(pncid,2)); %time. units: hours since 1800-01-01 00:00:0.0'
    ptime = (ptime/24)-21504; %dividing ptime by 24 gives the number of days since
    % 1800-01-01 00:00:0.0'. MJD is then obtained by subtracting 21504 (days between 
    % 01-01-1800 00:00:0.0' and 11-17-1858 00:00:0.0')
    pprecip = double(netcdf.getVar(pncid,3)); %precip data
    pprecip(pprecip < 0 ) = NaN;
    
    ptime_conv = datetime(ptime,'ConvertFrom','modifiedjuliandate');
    pyear_utc = ptime_conv.Year;
    pmonth_utc = ptime_conv.Month;
    pday_utc = ptime_conv.Day;
    phour_utc = (ptime_conv.Hour)*100;
    netcdf.close(pncid)
    clear pncid
    cd ../
    
    pday_calcount = find(pday_utc); % gives the day of the year (1-365/366)
    
  
    [~ , oiy] = size(tcyroi);
    TCP_daily_grid = zeros(length(plon),length(plat),length(pday_calcount)); %master TC precip grid
    [plonx , ~] = size(plon); %clear plony
    [platx , ~] = size(plat); %clear platy
    rplon = repmat(plon,1,120); % repeat longitude (300x1) 120x to make a 300x120 matrix
    rplat = repmat(plat,1,300)'; % repeat latitude (120x1) 300x to make a 120x300 matrix. Transpose for distance function.
    for a = 1:oiy % essentially the number of storms based on the length of the year % 9 is Katrina in 2005
        svl1 = find(tclon_oi(:,a) < 0,1,'first');
        evl1 = find(tclon_oi(:,a) < 0, 1,'last'); 
        svl2 = find(tclat_oi(:,a) > 0, 1, 'first'); 
        evl2 = find(tclat_oi(:,a) > 0, 1, 'last'); 
        x = tclon_oi(svl1:evl1,a)';
        y = tclat_oi(svl2:evl2,a)';
     
        t = [0,cumsum(abs(diff(x)) + abs(diff(y)))];
        nancheck = isnan(x);
        if all(nancheck == 0) == 1 % this is the logical true (1), so while 
            % isnan gives you a logical matrix, the all statement is asking 
            % if all of the t is cumulative chordal arc length between the
            % points. Basically, this large if statement is designed to do
            % spatial interpolation between synoptic time intervals.
            interruptions = 0;
            if length(x) == length(unique(t)) && length(y) == length(unique(t))
                npts=2*length(x); % number of points to generate along storm trajectory.
                % It's a function of the storm duration.
                % Set them up as linearly spaced along trajectory
                ti = linspace(0,t(end),npts);
                
                % estimate x,y coordinates along the trajectory using a
                % cubic spline
                xi = spline(t,x,ti); xi = xi';
                yi = spline(t,y,ti); yi = yi';
            else
                if length(x) ~= length(unique(t)) && length(y) == length(unique(t))
                    while length(x) ~= length(unique(t)) % this while loop should replace
                        % everything between the while and its parent if statement
                        x = x';
                        [~, ind] = unique(x(:,:),'rows','stable');
                        duplicate_ind = setdiff(1:size(x,1), ind);
                        ldupin = length(duplicate_ind);
                        x(duplicate_ind, :) = x(duplicate_ind,:)+(rand(ldupin,1)-0.5)/1e6;
                        x = x';
                        t = [0, cumsum(abs(diff(x))+abs(diff(y)))];
                    end
                    npts = 2*length(x);
                    ti = linspace(0,t(end),npts);
                    
                    xi = spline(t,x,ti); xi = xi';
                    yi = spline(t,y,ti); yi = yi';
                    
                elseif length(x) == length(unique(t)) && length(y) ~= length(unique(t))
                    while length(y) ~= length(unique(t)) % same thing as the previous while loop and if statement
                        y=y';
                        [~,ind] = unique(y(:,:),'rows','stable');
                        duplicate_ind = setdiff(1:size(y,1), ind);
                        ldupin = length(duplicate_ind);
                        y(duplicate_ind,:)=y(duplicate_ind,:)+(rand(ldupin,1)-0.5)/1e6;
                        y=y';
                        t = [0, cumsum(abs(diff(x))+abs(diff(y)))];
                    end
                    npts = 2*length(x);
                    ti = linspace(0,t(end),npts);
                    
                    xi = spline(t,x,ti); xi = xi';
                    yi = spline(t,y,ti); yi = yi'; 
                    
                elseif length(x) ~= length(unique(t)) && length(y) ~= length(unique(t))
                    while length(x) ~= length(unique(t)) || length(y) ~= length(unique(t)) % and same thing, but more complicated
                        if length(x) ~= length(unique(t)) && length(y) ~= length(unique(t))
                            x=x';
                            [~, indx] = unique(x(:,:),'rows','stable');
                            duplicate_indx = setdiff(1:size(x,1), indx);
                            ldupinx = length(duplicate_indx);
                            x(duplicate_indx,:)=x(duplicate_indx,:)+(rand(ldupinx,1)-0.5)/1e6;
                            x=x';
                            y=y';
                            [~, indy] = unique(y(:,:),'rows','stable');
                            duplicate_indy = setdiff(1:size(y,1), indy);
                            ldupiny = length(duplicate_indy);
                            y(duplicate_indy,:)=y(duplicate_indy,:)+(rand(ldupiny,1)-0.5)/1e6;
                            y=y';
                            t = [0, cumsum(abs(diff(x))+abs(diff(y)))];
                        elseif length(x) ~= length(unique(t)) && length(y) == length(unique(t))
                            x=x';
                            [~, indx] = unique(x(:,:),'rows','stable');
                            duplicate_indx = setdiff(1:size(x,1), indx);
                            ldupinx = length(duplicate_indx);
                            x(duplicate_indx, :)=x(duplicate_indx,:)+(rand(ldupinx,1)-0.5)/1e6;
                            x=x';
                            t = [0, cumsum(abs(diff(x))+abs(diff(y)))];
                        elseif length(x) == length(unique(t)) && length(y) ~= length(unique(t))
                            y=y';
                            [~, indy] = unique(y(:,:),'rows','stable');
                            duplicate_indy = setdiff(1:size(y,1), indy);
                            ldupiny = length(duplicate_indy);
                            y(duplicate_indy,:)=y(duplicate_indy,:)+(rand(ldupiny,1)-0.5)/1e6;
                            y=y';
                            t = [0, cumsum(abs(diff(x))+abs(diff(y)))];
                        end
                    end
                    npts = 2*length(x);
                    ti = linspace(0,t(end),npts);
                    
                    xi = spline(t,x,ti); xi = xi';
                    yi = spline(t,y,ti); yi = yi';
                end
            end
            
            
            hurhour_utc = repelem(tchour_oi(svl1:evl1,a),2); % TC UTC hour and double it
            hurhour_utc_loc = find(hurhour_utc >= 0);
            for hr = 1:length(hurhour_utc)
                if mod(hurhour_utc_loc(hr,1),2) == 0
                    hurhour_utc(hr,1) = hurhour_utc(hr,1)+300;
                else
                    hurhour_utc(hr,1) = hurhour_utc(hr,1);
                end
            end
            hur_doy = repelem(tcdate_doy(svl1:evl1,a),2); % doubling each DotY value to match the interpolated coordinates
            for hr = 1:length(hurhour_utc)
                if hurhour_utc(hr,1) > 1200 % Precip is accumulated from -12Z to +12Z.
                    hur_doy(hr,1) = hur_doy(hr,1)+1;
                else
                    hur_doy(hr,1) = hur_doy(hr,1);
                end
            end
            hdoy_uniq = unique(hur_doy); % unique DotY
            hdoy_CONUS_precip = pprecip(:,:,hdoy_uniq); % CONUS precipitation on the unique TC days
            % hdoy_grid_boo is the Boolean array saying whether a TC was there (1)
            % or not (0). It resets with each iteration of a because I have issues.
            hdoy_grid_boo = zeros(length(plon),length(plat),length(hdoy_uniq));
            [hdoyux , ~] = size(hdoy_uniq);
        elseif any(nancheck) == 0 % add the logical false (0) if there are any NaNs.
            interruptions = 1;
            prenan_loc = find(isnan(tclon_oi(:,a)),1,'first')-1;
            postnan_loc = find(isnan(tclon_oi(:,a)),1,'last')+1;
            x = tclon_oi(svl1:prenan_loc,a)';
            x2 = tclon_oi(postnan_loc:evl1,a)';
            y = tclat_oi(svl1:prenan_loc,a)';
            y2 = tclat_oi(postnan_loc:evl1,a)';
            t = [0, cumsum(abs(diff(x))+abs(diff(y)))];
            t2 = [0, cumsum(abs(diff(x2))+abs(diff(y2)))];
            if length(x) == length(unique(t)) && length(y) == length(unique(t))
                npts=2*length(x);
                
                % set them up as linearly spaced along trajectory
                ti = linspace(0,t(end),npts);
                
                % estimate x,y coordinates along the trajectory
                xi = spline(t,x,ti); xi = xi';
                yi = spline(t,y,ti); yi = yi';
                
            else
                if length(x) ~= length(unique(t)) && length(y) == length(unique(t))
                    while length(x) ~= length(unique(t)) % You get it. I'm repeating things.
                        x = x';
                        [~, ind] = unique(x(:,:),'rows','stable');
                        duplicate_ind = setdiff(1:size(x,1), ind);
                        ldupin = length(duplicate_ind);
                        x(duplicate_ind, :) = x(duplicate_ind,:)+(rand(ldupin,1)-0.5)/1e6;
                        x = x';
                        t = [0, cumsum(abs(diff(x))+abs(diff(y)))];
                    end
                    npts = 2*length(x);
                    ti = linspace(0,t(end),npts);
                    
                    xi = spline(t,x,ti); xi = xi';
                    yi = spline(t,y,ti); yi = yi';
                    
                elseif length(x) == length(unique(t)) && length(y) ~= length(unique(t))
                    while length(y) ~= length(unique(t)) % yep
                        y=y';
                        [~,ind] = unique(y(:,:),'rows','stable');
                        duplicate_ind = setdiff(1:size(y,1), ind);
                        ldupin = length(duplicate_ind);
                        y(duplicate_ind,:)=y(duplicate_ind,:)+(rand(ldupin,1)-0.5)/1e6;
                        y=y';
                        t = [0, cumsum(abs(diff(x))+abs(diff(y)))];
                    end

                    npts = 2*length(x);
                    ti = linspace(0,t(end),npts);
                    
                    xi = spline(t,x,ti); xi = xi';
                    yi = spline(t,y,ti); yi = yi';
               
                elseif length(x) ~= length(unique(t)) && length(y) ~= length(unique(t))
                    while length(x) ~= length(unique(t)) || length(y) ~= length(unique(t)) % same
                        if length(x) ~= length(unique(t)) && length(y) ~= length(unique(t))
                            x=x';
                            [~, indx] = unique(x(:,:),'rows','stable');
                            duplicate_indx = setdiff(1:size(x,1), indx);
                            ldupinx = length(duplicate_indx);
                            x(duplicate_indx,:)=x(duplicate_indx,:)+(rand(ldupinx,1)-0.5)/1e6;
                            x=x';
                            y=y';
                            [~, indy] = unique(y(:,:),'rows','stable');
                            duplicate_indy = setdiff(1:size(y,1), indy);
                            ldupiny = length(duplicate_indy);
                            y(duplicate_indy,:)=y(duplicate_indy,:)+(rand(ldupiny,1)-0.5)/1e6;
                            y=y';
                            t = [0, cumsum(abs(diff(x))+abs(diff(y)))];
                        elseif length(x) ~= length(unique(t)) && length(y) == length(unique(t))
                            x=x';
                            [~, indx] = unique(x(:,:),'rows','stable');
                            duplicate_indx = setdiff(1:size(x,1), indx);
                            ldupinx = length(duplicate_indx);
                            x(duplicate_indx, :)=x(duplicate_indx,:)+(rand(ldupinx,1)-0.5)/1e6;
                            x=x';
                            t = [0, cumsum(abs(diff(x))+abs(diff(y)))];
                        elseif length(x) == length(unique(t)) && length(y) ~= length(unique(t))
                            y=y';
                            [~, indy] = unique(y(:,:),'rows','stable');
                            duplicate_indy = setdiff(1:size(y,1), indy);
                            ldupiny = length(duplicate_indy);
                            y(duplicate_indy,:)=y(duplicate_indy,:)+(rand(ldupiny,1)-0.5)/1e6;
                            y=y';
                            t = [0, cumsum(abs(diff(x))+abs(diff(y)))];
                        end
                    end
                    npts = 2*length(x);
                    ti = linspace(0,t(end),npts);
                    
                    xi = spline(t,x,ti); xi = xi';
                    yi = spline(t,y,ti); yi = yi';
                end
            end
            if length(x2) == length(unique(t2)) && length(y2) == length(unique(t2))
                npts2=2*length(x2);
                
                % set them up as linearly spaced along trajectory
                ti2 = linspace(0,t2(end),npts2);
                
                % estimate x,y coordinates along the trajectory
                xi2 = spline(t2,x2,ti2); xi2 = xi2';
                yi2 = spline(t2,y2,ti2); yi2 = yi2';
               
            else
                if length(x2) ~= length(unique(t2)) && length(y2) == length(unique(t2))
                    while length(x2) ~= length(unique(t2)) % same
                        x2 = x2';
                        [~, ind2] = unique(x2(:,:),'rows','stable');
                        duplicate_ind2 = setdiff(1:size(x2,1), ind2);
                        ldupin2 = length(duplicate_ind2);
                        x2(duplicate_ind2, :) = x2(duplicate_ind2,:)+(rand(ldupin2,1)-0.5)/1e6;
                        x2 = x2';
                        t2 = [0, cumsum(abs(diff(x2))+abs(diff(y2)))];
                    end
                    npts2 = 2*length(x2);
                    ti2 = linspace(0,t2(end),npts2);
                    
                    xi2 = spline(t2,x2,ti2); xi2 = xi2';
                    yi2 = spline(t2,y2,ti2); yi2 = yi2';
                    
                elseif length(x2) == length(unique(t2)) && length(y2) ~= length(unique(t2))
                    while length(y2) ~= length(unique(t2)) % this is getting repetitive
                        y2=y2';
                        [~,ind2] = unique(y2(:,:),'rows','stable');
                        duplicate_ind2 = setdiff(1:size(y2,1), ind2);
                        ldupin2 = length(duplicate_ind2);
                        y2(duplicate_ind2,:)=y2(duplicate_ind2,:)+(rand(ldupin2,1)-0.5)/1e6;
                        y2=y2';
                        t2 = [0, cumsum(abs(diff(x2))+abs(diff(y2)))];
                    end
                    
                    npts2 = 2*length(x2);
                    ti2 = linspace(0,t2(end),npts2);
                    
                    xi2 = spline(t2,x2,ti2); xi2 = xi2';
                    yi2 = spline(t2,y2,ti2); yi2 = yi2';
                    
                elseif length(x2) ~= length(unique(t2)) && length(y2) ~= length(unique(t2))
                    while length(x2) ~= length(unique(t2)) || length(y2) ~= length(unique(t2)) % okay, good
                        if length(x2) ~= length(unique(t2)) && length(y2) ~= length(unique(t2))
                            x2=x2';
                            [~, indx2] = unique(x2(:,:),'rows','stable');
                            duplicate_indx2 = setdiff(1:size(x2,1), indx2);
                            ldupinx2 = length(duplicate_indx2);
                            x2(duplicate_indx2,:)=x2(duplicate_indx2,:)+(rand(ldupinx2,1)-0.5)/1e6;
                            x2=x2';
                            y2=y2';
                            [~, indy2] = unique(y2(:,:),'rows','stable');
                            duplicate_indy2 = setdiff(1:size(y2,1), indy2);
                            ldupiny2 = length(duplicate_indy2);
                            y2(duplicate_indy2,:)=y2(duplicate_indy2,:)+(rand(ldupiny2,1)-0.5)/1e6;
                            y2=y2';
                            t2 = [0, cumsum(abs(diff(x2))+abs(diff(y2)))];
                        elseif length(x2) ~= length(unique(t2)) && length(y2) == length(unique(t2))
                            x2=x2';
                            [~, indx2] = unique(x2(:,:),'rows','stable');
                            duplicate_indx2 = setdiff(1:size(x2,1), indx2);
                            ldupinx2 = length(duplicate_indx2);
                            x2(duplicate_indx2, :)=x2(duplicate_indx2,:)+(rand(ldupinx2,1)-0.5)/1e6;
                            x2=x2';
                            t2 = [0, cumsum(abs(diff(x2))+abs(diff(y2)))];
                        elseif length(x2) == length(unique(t2)) && length(y2) ~= length(unique(t2))
                            y2=y2';
                            [~, indy2] = unique(y2(:,:),'rows','stable');
                            duplicate_indy2 = setdiff(1:size(y2,1), indy2);
                            ldupiny2 = length(duplicate_indy2);
                            y2(duplicate_indy2,:)=y2(duplicate_indy2,:)+(rand(ldupiny2,1)-0.5)/1e6;
                            y2=y2';
                            t2 = [0, cumsum(abs(diff(x2))+abs(diff(y2)))];
                        end
                    end
                    npts2 = 2*length(x2);
                    ti2 = linspace(0,t2(end),npts2);
                    
                    xi2 = spline(t2,x2,ti2); xi2 = xi2';
                    yi2 = spline(t2,y2,ti2); yi2 = yi2';

                end
            end
            hurhour_utc = repelem(tchour_oi(svl1:prenan_loc,a),2); % TC UTC hour and double it
            hurhour_utc_loc = find(hurhour_utc >= 0);
            for hr = 1:length(hurhour_utc)
                if mod(hurhour_utc_loc(hr,1),2) == 0
                    hurhour_utc(hr,1) = hurhour_utc(hr,1)+300;
                else
                    hurhour_utc(hr,1) = hurhour_utc(hr,1);
                end
            end
            hur_doy = repelem(tcdate_doy(svl1:prenan_loc,a),2); % doubling each DotY 
            % value to match the interpolated coordinates
            for hr = 1:length(hurhour_utc)
                if hurhour_utc(hr,1) > 1200 % Precip accumulated from -12Z to +12Z.
                    hur_doy(hr,1) = hur_doy(hr,1)+1;
                else
                    hur_doy(hr,1) = hur_doy(hr,1);
                end
            end
            
            hdoy_uniq = unique(hur_doy); % unique DotY
            hdoy_CONUS_precip = pprecip(:,:,hdoy_uniq); % CONUS precipitation on the unique TC days
            % hdoy_grid_boo is the Boolean array saying whether a TC was there (1)
            % or not (0). It resets with each iteration of a because I have issues.
            hdoy_grid_boo = zeros(length(plon),length(plat),length(hdoy_uniq));
            [hdoyux , ~] = size(hdoy_uniq);
             
            hurhour_utc2 = repelem(tchour_oi(postnan_loc:evl1,a),2); % TC UTC hour and double it
            hurhour_utc_loc2 = find(hurhour_utc2 >= 0);
            for hr2 = 1:length(hurhour_utc2)
                if mod(hurhour_utc_loc2(hr2,1),2) == 0
                    hurhour_utc2(hr2,1) = hurhour_utc2(hr2,1)+300;
                else
                    hurhour_utc2(hr2,1) = hurhour_utc2(hr2,1);
                end
            end
            hur_doy2 = repelem(tcdate_doy(postnan_loc:evl1,a),2);
            for hr2 = 1:length(hurhour_utc2)
                if hurhour_utc2(hr2,1) > 1200 % Precip accumulated from -12Z to +12Z
                    hur_doy2(hr2,1) = hur_doy2(hr2,1)+1;
                else
                    hur_doy2(hr2,1) = hur_doy2(hr2,1);
                end
            end
            hdoy_uniq2 = unique(hur_doy2);
            hdoy_CONUS_precip2 = pprecip(:,:,hdoy_uniq2);
            hdoy_grid_boo2 = zeros(length(plon),length(plat),length(hdoy_uniq2));
            [hdoyux2 , ~] = size(hdoy_uniq2);
        end
        
        
        % This for loop basically gets all of the cells that are within the TC
        % search radius. It then populates hdoy_grid_boo with 1s and 0s to
        % indicate if a TC was present or absent, respectively.
        if interruptions == 0
            for k = 1:hdoyux
                doh = hur_doy((find(hur_doy == hdoy_uniq(k),1,'first')):(find(hur_doy == hdoy_uniq(k),1,'last')));
                doh_loc = (find(hur_doy == hdoy_uniq(k),1,'first'):find(hur_doy == hdoy_uniq(k),1,'last'));
                [~ , dohly] = size(doh_loc); % switch ~ with dohlx if you mess up
                tc_boo = zeros(plonx,platx,dohly);
                mm = 1;
                for m = doh_loc(1):doh_loc(dohly)
                    rxi = repmat(xi(m,1),300,120);
                    ryi = repmat(yi(m,1),300,120);
                    hdist = deg2km(distance(ryi,rxi,rplat,rplon),'earth'); %arclen = distance(yi(m),xi(m),plat(j),plon(i)); 
                    % the line with arclen in it is from an earlier iteration; I'm keeping it in because I'm scared.
                    tccells = hdist <= 223; % matrix of 1s and 0s indicating 
                    % the presence (1) or absence (0) of a TC at a cell
                    tc_boo(:,:,mm) = tccells; % basically, this should make 
                    % an array of presence/absence values for TCs
                    mm = mm+1;
                end
                tc_boosum = nansum(tc_boo,3); % sum up all of the presence/absence 
                % data for the day (one for each coordinate/measurement time)
                tc_boosum(tc_boosum > 1) = 1; % prevents double counting 
                % precipitation values for a grid point for the day
                hdoy_grid_boo(:,:,k) = tc_boosum; % puts tc_boosum into a 3d array 
                % of Boolean values for the entire TC.
            end
    
            % This for loop uses the extracted TCP and Boolean grid to add to the
            % master TCP array. To prevent the removal of precipitation from
            % another storm occurring on the same date and outside of the search
            % radius of the storm (causing the previous value to be multiplied by
            % 0), the Boolean cells containing 0 and cells with precipitation from
            % another storm in the master array are used to set conditional
            % statements. From there, we just change the cell in the Boolean grid
            % to 1 and then continue with the element-wise multiplication. I hope
            % it works... Addendum: it works.
            
            check_tcp_doys = TCP_daily_grid(:,:,hdoy_uniq);
            check_tcp = any(check_tcp_doys(:)>0);
            if a == 1 || check_tcp == 0
                TCP_daily_grid(:,:,hdoy_uniq(:)) = hdoy_CONUS_precip(:,:,:).*hdoy_grid_boo(:,:,:);
            else
                for k = 1:hdoyux
                    check_date = TCP_daily_grid(:,:,hdoy_uniq(k));
                    check_date_boo = check_date > 0;
                    hdoy_grid_boo_sum = hdoy_grid_boo(:,:,k) + check_date_boo;
                    hdoy_grid_boo_sum(hdoy_grid_boo_sum > 1) = 1;
                    hdoy_grid_boo(:,:,k) = hdoy_grid_boo_sum;
                end
                TCP_daily_grid(:,:,hdoy_uniq(:)) = hdoy_CONUS_precip(:,:,:).*hdoy_grid_boo(:,:,:);
            end     
            
        else
            for k = 1:hdoyux
                doh = hur_doy((find(hur_doy == hdoy_uniq(k),1,'first')):(find(hur_doy == hdoy_uniq(k),1,'last')));
                doh_loc = (find(hur_doy == hdoy_uniq(k),1,'first'):find(hur_doy == hdoy_uniq(k),1,'last'));
                [~ , dohly] = size(doh_loc); %switch ~ with dohlx if you mess up
                tc_boo = zeros(plonx,platx,dohly);
                mm = 1;
                for m = doh_loc(1):doh_loc(dohly)
                    rxi = repmat(xi(m,1),300,120);
                    ryi = repmat(yi(m,1),300,120);
                    hdist = deg2km(distance(ryi,rxi,rplat,rplon),'earth'); %arclen = distance(yi(m),xi(m),plat(j),plon(i));
                    tccells = hdist <= 223; % matrix of 1s and 0s indicating the 
                    % presence (1) or absence (0) of a TC at a cell
                    tc_boo(:,:,mm) = tccells; % basically, this should make an 
                    % array of presence/absence values for TCs
                    mm = mm+1;
                end
                tc_boosum = nansum(tc_boo,3); % sum up all of the presence/absence 
                % data for the day (one for each coordinate/measurement time)
                tc_boosum(tc_boosum > 1) = 1; % prevents double counting 
                % precipitation values for a grid point for the day
                hdoy_grid_boo(:,:,k) = tc_boosum; % puts tc_boosum into a 
                % 3d array of Boolean values for the entire TC.
            end
        
            for k2 = 1:hdoyux2
                doh2 = hur_doy2((find(hur_doy2 == hdoy_uniq2(k2),1,'first')):(find(hur_doy2 == hdoy_uniq2(k2),1,'last')));
                doh_loc2 = (find(hur_doy2 == hdoy_uniq2(k2),1,'first'):find(hur_doy2 == hdoy_uniq2(k2),1,'last'));
                [~ , dohly2] = size(doh_loc2); % switch ~ with dohlx2 if you mess up
                tc_boo2 = zeros(plonx,platx,dohly2);
                mm2 = 1;
                for m2 = doh_loc2(1):doh_loc2(dohly2)
                    rxi2 = repmat(xi2(m2,1),300,120);
                    ryi2 = repmat(yi2(m2,1),300,120);
                    hdist2 = deg2km(distance(ryi2,rxi2,rplat,rplon),'earth'); %arclen = distance(yi(m),xi(m),plat(j),plon(i));
                    tccells2 = hdist2 <= 223; % matrix of 1s and 0s indicating 
                    % the presence (1) or absence (0) of a TC at a cell
                    tc_boo2(:,:,mm2) = tccells2; % basically, this should make 
                    % an array of presence/absence values for TCs
                    mm2 = mm2+1;
                end
                tc_boosum2 = nansum(tc_boo2,3); % sum up all of the 
                % presence/absence data for the day (one for each 
                % coordinate/measurement time)
                tc_boosum2(tc_boosum2 > 1) = 1; % prevents double 
                % counting precipitation values for a grid point 
                % for the day
                hdoy_grid_boo2(:,:,k) = tc_boosum2; % puts tc_boosum 
                % into a 3d array of Boolean values for the entire TC.
            end
            
            % This for loop uses the extracted TCP and Boolean grid to add to the
            % master TCP array. To prevent the removal of precipitation from
            % another storm occurring on the same date and outside of the search
            % radius of the storm (causing the previous value to be multiplied by
            % 0), the Boolean cells containing 0 and cells with precipitation from
            % another storm in the master array are used to set conditional
            % statements. From there, we just change the cell in the Boolean grid
            % to 1 and then continue with the element-wise multiplication. I hope
            % it works... Addendum: it works.
            check_tcp_doys = TCP_daily_grid(:,:,hdoy_uniq);
            check_tcp = any(check_tcp_doys(:)>0);
            if a == 1 || check_tcp == 0
                TCP_daily_grid(:,:,hdoy_uniq(:)) = hdoy_CONUS_precip(:,:,:).*hdoy_grid_boo(:,:,:);
            else
                for k = 1:hdoyux
                    check_date = TCP_daily_grid(:,:,hdoy_uniq(k));
                    check_date_boo = check_date > 0;
                    hdoy_grid_boo_sum = hdoy_grid_boo(:,:,k) + check_date_boo;
                    hdoy_grid_boo_sum(hdoy_grid_boo_sum > 1) = 1;
                    hdoy_grid_boo(:,:,k) = hdoy_grid_boo_sum;
                end
                TCP_daily_grid(:,:,hdoy_uniq(:)) = hdoy_CONUS_precip(:,:,:).*hdoy_grid_boo(:,:,:);
            end     
            check_tcp_doys = TCP_daily_grid(:,:,hdoy_uniq2);
            check_tcp = any(check_tcp_doys(:)>0);
            if a == 1 || check_tcp == 0
                TCP_daily_grid(:,:,hdoy_uniq2(:)) = hdoy_CONUS_precip2(:,:,:).*hdoy_grid_boo2(:,:,:);
            else
                for kk = 1:hdoyux2
                    check_date = TCP_daily_grid(:,:,hdoy_uniq2(kk));
                    check_date_boo = check_date > 0;
                    hdoy_grid_boo_sum = hdoy_grid_boo2(:,:,kk) + check_date_boo;
                    hdoy_grid_boo_sum(hdoy_grid_boo_sum > 1) = 1;
                    hdoy_grid_boo2(:,:,kk) = hdoy_grid_boo_sum;
                end
                TCP_daily_grid(:,:,hdoy_uniq2(:)) = hdoy_CONUS_precip2(:,:,:).*hdoy_grid_boo2(:,:,:);
            end
        end
    end
    cd TCP_data\
    
    
    fn = ['TCP_daily_',num2str(b+1947,4),'.mat'];
    save(fn, 'TCP_daily_grid', '-mat')
  
    ncid_tcp_daily = netcdf.create(['TCP_daily_',num2str(b+(st_yr-1),4),'.nc'],'NC_WRITE');
    dimidlat_daily = netcdf.defDim(ncid_tcp_daily,'latitude',120);
    dimidlon_daily = netcdf.defDim(ncid_tcp_daily,'longitude',300);
    dimidt_daily = netcdf.defDim(ncid_tcp_daily,'time',length(pday_calcount));
    
    latitude_ID_daily = netcdf.defVar(ncid_tcp_daily,'latitude','NC_FLOAT',[dimidlat_daily]);
    longitude_ID_daily = netcdf.defVar(ncid_tcp_daily,'longitude','NC_FLOAT',[dimidlon_daily]);
    date_ID_daily = netcdf.defVar(ncid_tcp_daily,'time','NC_DOUBLE',[dimidt_daily]);
    TCPDat_daily_ID = netcdf.defVar(ncid_tcp_daily,'TCP_daily_grid','NC_FLOAT',[dimidlon_daily,dimidlat_daily,dimidt_daily]);
    netcdf.endDef(ncid_tcp_daily);
    netcdf.putVar(ncid_tcp_daily,longitude_ID_daily,plon);
    netcdf.putVar(ncid_tcp_daily,latitude_ID_daily,plat);
    netcdf.putVar(ncid_tcp_daily,date_ID_daily,pday_calcount);
    netcdf.putVar(ncid_tcp_daily,TCPDat_daily_ID,TCP_daily_grid);
    netcdf.close(ncid_tcp_daily);
    cd ../../
    TCP_ann_sum = sum(TCP_daily_grid,3);
    TCP_ann_sum(TCP_ann_sum == 0) = NaN;
    TCP_annual(:,:,b) = TCP_ann_sum;
    
end
cd TCP_data\
save TCP_annual.mat TCP_annual -mat
% Add code to create netcdf files.

ncid_tcp_ann = netcdf.create('TCPDat_annual.nc','NC_WRITE'); %Open the file
% Define dimensions
dimidlat_ann = netcdf.defDim(ncid_tcp_ann,'latitude',120); % This defines latitude. It should always be 120.
dimidlon_ann = netcdf.defDim(ncid_tcp_ann,'longitude',300); % This defines longitude. It should always be 300.
dimidt_ann = netcdf.defDim(ncid_tcp_ann,'time',length(yr_uniq)); % This is the date dimension. It should be based on the number of years the user is extracting.
% Define the IDs for dimension variable
latitude_ID_ann = netcdf.defVar(ncid_tcp_ann,'latitude','NC_FLOAT',[dimidlat_ann]); % Single datatype is used for lat and lon just because that's what URD uses.
longitude_ID_ann = netcdf.defVar(ncid_tcp_ann,'longitude','NC_FLOAT',[dimidlon_ann]);
date_ID_ann = netcdf.defVar(ncid_tcp_ann,'time','NC_DOUBLE',[dimidt_ann]); % Double is used because that's what URD uses
% Define main variable
TCPDat_annual_ID = netcdf.defVar(ncid_tcp_ann,'TCP_annual','NC_FLOAT',[dimidlon_ann,dimidlat_ann,dimidt_ann]); % Single is used like URD. Also, the dimensions are lon x lat x time, or 300x120x(user-defined time).
% Stop defining variables
netcdf.endDef(ncid_tcp_ann);
% Store dimension variables
netcdf.putVar(ncid_tcp_ann,longitude_ID_ann,plon);
netcdf.putVar(ncid_tcp_ann,latitude_ID_ann,plat);
year_netcdf_range = [st_yr:1:ed_yr]';
netcdf.putVar(ncid_tcp_ann,date_ID_ann,year_netcdf_range);
% Store the main variable
netcdf.putVar(ncid_tcp_ann,TCPDat_annual_ID,TCP_annual);
%close netcdf
netcdf.close(ncid_tcp_ann);
cd ../
