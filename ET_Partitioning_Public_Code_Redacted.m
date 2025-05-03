%% ET Partitioning Code 
%  Code for isotope-based ET partitioning analysis at 13 NEON sites. 
%  Manuscript Title: Evapotranspiration Partitioning Across US Ecoregions: a Multi-Site Study Using Field Stable-Isotope Observations
%  Date last updated: 05/03/2025

%% Gathering Data 
%Path to directory where the code is stored
%Example: 'C:\...\Code';
code_path = '';
cd (code_path)

%Path to directory where data products are is stored
%Example: 'C:\...\Data';
data_path = '';

%Load NEON sites Metadata Including ID, Sample Date, and NEON Domains 
siteid      = readtable("siteid.xlsx");
SampDate    = readtable("SampDate_start_and_end_day.xlsx");
site_date   = readtable("site_date.xlsx");
site_domain = readmatrix("site_domain.xlsx");

%% Variable: Rs       
% Description: Used to determine sunrise and sunset times (R>10W/m2)

% NEON Data Product used below: Shortwave radiation (direct and diffuse pyranometer; DP1.00014.001)
    % Set path to location of this NEON Shortwave radiation data product on the local computer
    NEON_Rs_path = 'NEON_rad-short-direct-diffuse';
    fullfile(data_path,NEON_Rs_path)
    cd (fullfile(data_path,NEON_Rs_path))

%NEON file naming convention
filname1    = 'NEON.D%02.f.%s.DP1.00014.001.%d-%02d*'; filname1x    = '*_30min.*';
pat_datestr = digitsPattern(4) + "-" + digitsPattern(2) + "-" + digitsPattern(2); 
pat_time    = digitsPattern(2) + ":" + digitsPattern(2)  + ":" + digitsPattern(2); 

%Adjusted dates due to noisy data
Samp_adjs = zeros(57,1); %empty matrix
Samp_adjs (7) = 1; Samp_adjs (11) = 11; Samp_adjs (13) = 16; Samp_adjs (19) = 10; Samp_adjs (21) = 1; Samp_adjs (33) = -4; Samp_adjs (45) = -7; Samp_adjs (53) = 7; Samp_adjs (54) = 4;

%Loop through sample dates, extract correct NEON_rad-short-direct-diffuse file
for u = 1 :size(SampDate,1)

    clear filname2 filname filnamex1 data SS datestrSS timeeSS data datax datax_sub idx daytime_Rs
    filname2    = sprintf (filname1, site_domain(u), siteid.site_name{u},site_date.Year(u), site_date.Month(u));
    filname     = dir(fullfile(pwd,filname2));
    
    cd(filname.name);
    filnamex1   =  dir(fullfile(pwd,filname1x));
    data        = readtable(filnamex1.name);
    
    SS          = data.startDateTime; 
    
    datestrSS   = extract(SS,pat_datestr);
    timeeSS     = extract(SS,pat_time);
    testSS      = datetime(datestrSS + " " + timeeSS, "InputFormat","yyyy-MM-dd HH:mm:ss","Format","dd-MMM-yyyy HH:mm:ss.SS");
    datax       = table2timetable(data,'RowTimes',testSS);
    
    %Sampdate in this case is the sample date (UTC) entire day
    
        if ismember(u,[7 11 13 19 21 33 45 53 54])
        
            % SampDate_start_end has the 24-hr samp period (UTC)
            datax_sub  = datax(timerange(SampDate.Start_Sample_Day(u)+ days(Samp_adjs(u)),SampDate.End_Sample_Day(u)+ days(Samp_adjs(u)),'closed'),:);
            idx        = find(datax_sub.dirRadMean > 10);
            daytime_Rs = datax_sub(idx,:);
        
            %SampDate_sunrise_set has the time (UTC) end and start of Rs >10 W/m2
            SampDate_sunrise_set (u,1) = daytime_Rs.Time(1) - days(Samp_adjs(u));
            SampDate_sunrise_set (u,2) = daytime_Rs.Time(end) - days(Samp_adjs(u));
        
        else 
            % SampDate_start_end has the 24-hr samp period (UTC)
            datax_sub  = datax(timerange(SampDate.Start_Sample_Day(u)-3,SampDate.End_Sample_Day(u)+3,'closed'),:);
            idx        = find(datax_sub.dirRadMean > 10);
            daytime_Rs = datax_sub(idx,:);
            
            %SampDate_sunrise_set has the time (UTC) end and start of Rs >10 W/m2
            SampDate_sunrise_set (u,1) = daytime_Rs.Time(1);
            SampDate_sunrise_set (u,2) = daytime_Rs.Time(end);
        
        end
    cd(fullfile(code_path))
    writetimetable(daytime_Rs,'Rs_Matlab.xlsx','Sheet',siteid.site_id{u});
    cd (fullfile(data_path,NEON_Rs_path))
   

end

%Save output
cd(fullfile(code_path))
writematrix(SampDate_sunrise_set,'Rs_Matlab.xlsx','Sheet','Rs_Summarized');
writematrix(SampDate_sunrise_set,'SampDate_sunrise_set_Matlab.xlsx');

%% Variable: Ts       
% Description: Soil temperature

cd (code_path)

%Individual sensor naming convention
Ts_string = readtable('Ts_string.xlsx');

% NEON Data Product used below: Soil temperature (DP1.00041.001).
    % Set path to location of this NEON Soil Temperature data product on the local computer
    NEON_Ts_path = 'NEON_temp-soil';
    cd (fullfile(data_path,NEON_Ts_path))

%NEON file naming conventino
filname1    = 'NEON.D%02.f.%s.DP1.00041.001.%d-%02d*';
filname1x    = '*.501.030.*'; filname2x = '*.502.030.*'; filname3x = '*.503.030.*';
pat_datestr = digitsPattern(4) + "-" + digitsPattern(2) + "-" + digitsPattern(2); 
pat_time    = digitsPattern(2) + ":" + digitsPattern(2)  + ":" + digitsPattern(2); 

%Loop through sample dates, extract correct NEON Soil Temperature data file
for u = 1:length(siteid.site_id)

filname2    = sprintf (filname1, site_domain(u), siteid.site_name{u},site_date.Year(u), site_date.Month(u));
filname     = dir(fullfile(pwd,filname2));

cd(filname.name)

filnamex1    = struct2table(dir(fullfile(pwd,filname1x)));
filnamex2    = struct2table(dir(fullfile(pwd,filname2x)));
filnamex3    = struct2table(dir(fullfile(pwd,filname3x)));
filnamex     = vertcat(filnamex1,filnamex2,filnamex3);

clear currdepth currdepth_TT
dd = 0;
currdepth_str = cell(1,15);

%Loop through each sensor for any given site/sample date
%currdepth stores the data for  given sensor for the sample date of interest
%currdepth_str stores the sensor ID/name
for z = 1:size(filnamex,1)
    clear data SS datestrSS timeeSS datax datax_sub
    data        = readtable(filnamex.name{z});
    SS          = data.startDateTime; 
    datestrSS   = extract(SS,pat_datestr);
    timeeSS     = extract(SS,pat_time);
    testSS      = datetime(datestrSS + " " + timeeSS, "InputFormat","yyyy-MM-dd HH:mm:ss","Format","dd-MMM-yyyy HH:mm:ss.SS");
    datax       = table2timetable(data,'RowTimes',testSS);
    datax_sub   = datax(timerange(SampDate_sunrise_set(u,1),SampDate_sunrise_set(u,2),'closed'),:);
    if ~isnan(sum(datax_sub.soilTempMean)) 
        dd = dd+1;
        currdepth(:,dd) = datax_sub.soilTempMean;
        currdepth_str{1,dd} =  convertCharsToStrings(Ts_string.str{dd});
    end
end

currdepth_TT   = array2timetable(currdepth,'RowTimes',datax_sub.Time);

%Calculate and save the mean, std. dev, and number of sensors used
Ts_summ_out(u,1) = mean(mean(currdepth,1,'omitnan'));
Ts_summ_out(u,2) = std(mean(currdepth,1,'omitnan'))/sqrt(size(currdepth,2));
Ts_summ_out(u,3) = size(currdepth,2);

cd(fullfile(code_path))
writetimetable(currdepth_TT,'Ts_Matlab.xlsx','Sheet',siteid.site_id{u});
writecell(currdepth_str,'Ts_Matlab_CorrSP_and_Depths_Used.xlsx','Sheet',siteid.site_id{u});
cd (fullfile(data_path,NEON_Ts_path))

end

%Save output
cd(fullfile(code_path))
writematrix(Ts_summ_out,'Ts_Matlab.xlsx','Sheet','Ts_Summarized');


%% Variable: Ta       
% Description: Air Temperature

cd (code_path)

%Number of levels at each Eddy Covriance tower
eddy_levels = readmatrix('eddy_levels.xlsx');

% NEON Data Product used below: Bundled data products - eddy covariance (DP4.00200.001).
    % Set path to location of this NEON data product on the local computer
    NEON_Eddy_path = 'NEON-eddy-covariance';
    cd (fullfile(data_path,NEON_Eddy_path))

%NEON file naming conventino
filname1 = 'NEON.D%02.f.%s.DP4.00200.001.%d-%02d*';
pat      = digitsPattern(4) + "-" + digitsPattern(2);

%Loop through sample dates, extract correct NEON Air Temperature data file
for u = 1:length(siteid.site_id)

clear t site filname2 filname h5filezipped h5temp info_h5 ext_temp yr_temp month_temp
site       = string(siteid.site_name{u});
filname2   = sprintf (filname1, site_domain(u), siteid.site_name{u},site_date.Year(u), site_date.Month(u));
filname    = dir(fullfile(pwd,filname2));

cd(filname.name)
h5temp = dir ('*.h5');
info_h5 = h5info(h5temp.name);

%Extract date from file since time begin and time end is unreadable in the h5df file
ext_temp   = extract(string(filname.name),pat);
dummy_347  = split(ext_temp,'-');
yr_temp    = str2double(dummy_347(1));
month_temp = str2double(dummy_347(2));
t1         = datetime(yr_temp,month_temp,1);
tend       = datetime(yr_temp,month_temp,eomday(yr_temp,month_temp),23,30,0);
tstep      = minutes(30);
t          = (t1:tstep:tend)';

%Loop through levels of Eddy Tower and extract air temperature
clear data_tempAirLvl datax datax_sub currLvl currLvl_TT 
dd = 0; %dummy variable

for l = 2:2:eddy_levels(u)*2 
    clear path_tempAirLvl path_tempAirLvl_temp data_tempAirLvl datax datax_sub 
if (l/2 ==   eddy_levels(u)) == 0
    path_tempAirLvl      = getfield(info_h5(1).Groups.Groups(1).Groups(1).Groups(14).Groups,{l},'Name');
    path_tempAirLvl_temp = append(path_tempAirLvl,'/temp');
    data_tempAirLvl      = h5read(h5temp.name,path_tempAirLvl_temp);

    datax       = array2timetable(data_tempAirLvl.mean,'RowTimes',t);
    datax_sub   = datax(timerange(SampDate.Start_Sample_Day(u),SampDate.End_Sample_Day(u),'closed'),:);

if any(sum(datax_sub.Var1,"omitnan"))  
    dd = dd+1;
    currLvl(:,dd) = datax_sub.Var1;
    continue
end

else
    path_tempAirLvl      = getfield(info_h5(1).Groups.Groups(1).Groups(1).Groups(15).Groups,{2},'Name');
    path_tempAirLvl_temp = append(path_tempAirLvl,'/temp');
    data_tempAirLvl      = h5read(h5temp.name,path_tempAirLvl_temp);

    datax       = array2timetable(data_tempAirLvl.mean,'RowTimes',t);
    datax_sub   = datax(timerange(SampDate.Start_Sample_Day(u),SampDate.End_Sample_Day(u),'closed'),:);

if any(sum(datax_sub.Var1,"omitnan")) 
    dd = dd+1;
    currLvl(:,dd) = datax_sub.Var1;
    continue
end
end
end

%Calculate and save the mean, std. dev, and number of Eddy Covariance levels used
currLvl_TT   = array2timetable(currLvl,'RowTimes',datax_sub.Time);
currLvl_summ_out(u,1) = mean(mean(currLvl,1,'omitnan'),'omitnan');
currLvl_summ_out(u,2) = std(mean(currLvl,1,'omitnan'),'omitnan')/sqrt(size(currLvl,2));
currLvl_summ_out(u,3) = size(currLvl,2);

cd(fullfile(code_path))
writetimetable(currLvl_TT,'Ta_Matlab.xlsx','Sheet',siteid.site_id{u});
cd (fullfile(data_path,NEON_Eddy_path))

end

%Save output
cd(fullfile(code_path))
writematrix(currLvl_summ_out,'Ta_Matlab.xlsx','Sheet','Ta_Summarized');
writematrix(currLvl_summ_out,'Ta_summ_out_Matlab.xlsx');

%% Variable: RH       
% Description: Relative Humidity

% NEON Data Product used below: Bundled data products - eddy covariance (DP4.00200.001).
    % Set path to location of this NEON data product on the local computer
    NEON_Eddy_path = 'NEON-eddy-covariance';
    cd (fullfile(data_path,NEON_Eddy_path))

%NEON file naming conventino
filname1    = 'NEON.D%02.f.%s.DP4.00200.001.%d-%02d*';
pat        = digitsPattern(4) + "-" + digitsPattern(2);

skipu = [13 20 21 37];

%Loop through sample dates, extract correct NEON data file
for u = 1:length(siteid.site_id)

    if ~ismember(u,skipu) 

site = string(siteid.site_name{u});
filname2    = sprintf (filname1, site_domain(u), siteid.site_name{u},site_date.Year(u), site_date.Month(u));
filname     = dir(fullfile(pwd,filname2));

cd(filname.name)
h5temp = dir ('*.h5');
info_h5 = h5info(h5temp.name);

%Extract date from file since time begin and time end is nonsense in the h5df file
ext_temp   = extract(string(filname.name),pat);
dummy_347  = split(ext_temp,'-');
yr_temp    = str2double(dummy_347(1));
month_temp = str2double(dummy_347(2));
t1         = datetime(yr_temp,month_temp,1);
tend       = datetime(yr_temp,month_temp,eomday(yr_temp,month_temp),23,30,0);
tstep      = minutes(30);
t          = (t1:tstep:tend)';

clear data_tempAirLvl datax datax_sub currLvl currLvl_TT 
dd = 0;
for l = 2:2:eddy_levels(u)*2 %Loop through levels of Eddy Towers 
    clear path_tempAirLvl path_tempAirLvl_temp data_tempAirLvl datax datax_sub 
    path_tempRHLvl      = getfield(info_h5(1).Groups.Groups(1).Groups(1).Groups(10).Groups,{l},'Name');
    path_tempRHLvl_temp = append(path_tempRHLvl,'/rhEnvHut');
    data_tempRHLvl      = h5read(h5temp.name,path_tempRHLvl_temp);

    datax       = array2timetable(data_tempRHLvl.mean,'RowTimes',t);
    datax_sub   = datax(timerange(SampDate.Start_Sample_Day(u),SampDate.End_Sample_Day(u),'closed'),:);

if any(sum(datax_sub.Var1,"omitnan"))  
    dd = dd+1;
    currLvl(:,dd) = datax_sub.Var1;
    continue
end

end

%Calculate and save the mean, std. dev, and number of Eddy Covariance levels used
currLvl_TT   = array2timetable(currLvl,'RowTimes',datax_sub.Time);
RH_summ_out(u,1) = mean(mean(currLvl,1,'omitnan'),'omitnan');
RH_summ_out(u,2) = std(mean(currLvl,1,'omitnan'),'omitnan')/sqrt(size(currLvl,2));
RH_summ_out(u,3) = size(currLvl,2);

cd(fullfile(code_path))
writetimetable(currLvl_TT,'RH_Matlab.xlsx','Sheet',siteid.site_id{u});
cd (fullfile(data_path,NEON_Eddy_path))

    end

end

%Save output
cd(fullfile(code_path))
writematrix(RH_summ_out,'RH_Matlab.xlsx','Sheet','Ta_Summarized');
writematrix(RH_summ_out,'RH_summ_out_Matlab.xlsx');

%% Variable: dA (18O) 
% Description: Isotopic composition of atmospheric water vapor (18O)

% NEON Data Product used below: Bundled data products - eddy covariance (DP4.00200.001).
    % Set path to location of this NEON data product on the local computer
    NEON_Eddy_path = 'NEON-eddy-covariance';
    cd (fullfile(data_path,NEON_Eddy_path))

%NEON file naming conventino
filname1    = 'NEON.D%02.f.%s.DP4.00200.001.%d-%02d*';
pat        = digitsPattern(4) + "-" + digitsPattern(2);

%Loop through sample dates, extract correct NEON data file
for u = 1:length(siteid.site_id)

clear t site filname2 filname h5filezipped h5temp info_h5 ext_temp yr_temp month_temp currLvl_TT
site = string(siteid.site_name{u});
filname2    = sprintf (filname1, site_domain(u), siteid.site_name{u},site_date.Year(u), site_date.Month(u));
filname     = dir(fullfile(pwd,filname2));

cd(filname.name)

h5temp = dir ('*.h5');
info_h5 = h5info(h5temp.name);

%Extract date from file since time begin and time end is nonsense in the h5df file
ext_temp   = extract(string(filname.name),pat);
dummy_347  = split(ext_temp,'-');
yr_temp    = str2double(dummy_347(1));
month_temp = str2double(dummy_347(2));
t1         = datetime(yr_temp,month_temp,1);
tend       = datetime(yr_temp,month_temp,eomday(yr_temp,month_temp),23,30,0);
tstep      = minutes(30);
t          = (t1:tstep:tend)';

clear data_tempRHLvl datax datax_sub currLvl currLvl_TT
dd=0;
for l = 2:2:eddy_levels(u)*2 %Loop through levels of Eddy Towers 

clear path_tempRHLvl path_tempdALvl_18O_temp data_tempRHLvl datax datax_sub 

path_tempdALvl      = getfield(info_h5(1).Groups.Groups(1).Groups(1).Groups(10).Groups,{l},'Name');
path_tempdALvl_18O_temp = append(path_tempdALvl,'/dlta18OH2o');

data_tempdA_18O_Lvl      = h5read(h5temp.name,path_tempdALvl_18O_temp);

datax       = array2timetable(data_tempdA_18O_Lvl.mean,'RowTimes',t);
datax_sub   = datax(timerange(SampDate.Start_Sample_Day(u),SampDate.End_Sample_Day(u),'closed'),:);

if any(sum(datax_sub.Var1,"omitnan"))  
    dd = dd+1;
    currLvl(:,dd) = datax_sub.Var1;
    continue
end

end

currLvl_TT   = array2timetable(currLvl,'RowTimes',datax_sub.Time);

%Calculate and save the mean, std. dev, and number of Eddy Covariance levels used
dA_18O_Lvl_summ_out(u,1) = mean(mean(currLvl,1,'omitnan'),'omitnan');
dA_18O_Lvl_summ_out(u,2) = std(mean(currLvl,1,'omitnan'),'omitnan')/sqrt(size(currLvl,2));
dA_18O_Lvl_summ_out(u,3) = size(currLvl,2);

cd(fullfile(code_path))
writetimetable(currLvl_TT,'dA_18O_Matlab.xlsx','Sheet',siteid.site_id{u});
cd (fullfile(data_path,NEON_Eddy_path))

end

%Save output
cd(fullfile(code_path))
writematrix(dA_18O_Lvl_summ_out,'dA_18O_Matlab.xlsx','Sheet','dA_18O_Summarized');
writematrix(dA_18O_Lvl_summ_out,'dA_18O_Lvl_summ_out.xlsx','Sheet','dA_18O_Summarized');

%% Variable: dA (2H)  
% Description: Isotopic composition of atmospheric water vapor (2H)

% NEON Data Product used below: Bundled data products - eddy covariance (DP4.00200.001).
    % Set path to location of this NEON data product on the local computer
    NEON_Eddy_path = 'NEON-eddy-covariance';
    cd (fullfile(data_path,NEON_Eddy_path))

%NEON file naming convention
filname1    = 'NEON.D%02.f.%s.DP4.00200.001.%d-%02d*';
pat        = digitsPattern(4) + "-" + digitsPattern(2);

%Loop through sample dates, extract correct NEON data file
for u = 1:length(siteid.site_id)

clear t site filname2 filname h5filezipped h5temp info_h5 ext_temp yr_temp month_temp currLvl_TT
site = string(siteid.site_name{u});
filname2    = sprintf (filname1, site_domain(u), siteid.site_name{u},site_date.Year(u), site_date.Month(u));
filname     = dir(fullfile(pwd,filname2));

cd(filname.name)
h5temp = dir ('*.h5');
info_h5 = h5info(h5temp.name);

%Extract date from file since time begin and time end is nonsense in the h5df file
ext_temp   = extract(string(filname.name),pat);
dummy_347      = split(ext_temp,'-');
yr_temp    = str2double(dummy_347(1));
month_temp = str2double(dummy_347(2));
t1      = datetime(yr_temp,month_temp,1);
tend    = datetime(yr_temp,month_temp,eomday(yr_temp,month_temp),23,30,0);
tstep       = minutes(30);
t = (t1:tstep:tend)';

clear data_tempRHLvl datax datax_sub currLvl currLvl_TT
dd=0;
for l = 2:2:eddy_levels(u)*2 %Loop through levels of Eddy Towers 

clear path_tempRHLvl path_tempdALvl_2H_temp data_tempRHLvl datax datax_sub 

path_tempdALvl      = getfield(info_h5(1).Groups.Groups(1).Groups(1).Groups(10).Groups,{l},'Name');
path_tempdALvl_2H_temp = append(path_tempdALvl,'/dlta2HH2o');

data_tempdA_2H_Lvl      = h5read(h5temp.name,path_tempdALvl_2H_temp);

datax       = array2timetable(data_tempdA_2H_Lvl.mean,'RowTimes',t);
datax_sub   = datax(timerange(SampDate.Start_Sample_Day(u),SampDate.End_Sample_Day(u),'closed'),:);

if any(sum(datax_sub.Var1,"omitnan"))  
    dd = dd+1;
    currLvl(:,dd) = datax_sub.Var1;
    continue
end

end

currLvl_TT   = array2timetable(currLvl,'RowTimes',datax_sub.Time);

%Calculate and save the mean, std. dev, and number of Eddy Covariance levels used
dA_2H_Lvl_summ_out(u,1) = mean(mean(currLvl,1,'omitnan'),'omitnan');
dA_2H_Lvl_summ_out(u,2) = std(mean(currLvl,1,'omitnan'),'omitnan')/sqrt(size(currLvl,2));
dA_2H_Lvl_summ_out(u,3) = size(currLvl,2);

cd(fullfile(code_path))
writetimetable(currLvl_TT,'dA_2H_Matlab.xlsx','Sheet',siteid.site_id{u});
cd (fullfile(data_path,NEON_Eddy_path))

end

%Save output
cd(fullfile(code_path))
writematrix(dA_2H_Lvl_summ_out,'dA_2H_Matlab.xlsx','Sheet','dA_2H_Summarized');
writematrix(dA_2H_Lvl_summ_out,'dA_2H_Lvl_summ_out.xlsx','Sheet','dA_2H_Summarized');

%% Variable: dT and dS
% Description: Isotopic composition of plant-water and soil-water

% Data Product used below: WaterIsotope Database, Project ID: 00384. 
    % Set path to location of this data product on the local computer
    %Example: 'C:\\WB';
    WaterIso_path = '';
    cd (fullfile(WaterIso_path))

WIDB_Table              = readtable('.csv'); %Water Isotopes Database CSV

%Loop through sample dates, extract correct site/sample date
for u = 1:length(siteid.site_id)

clear matching WIDB_SiteSub_Table WIDB_SiteSub_Soil_Table WIDB_SiteSub_Stem_Table WIDB_SiteSub_dSoil_TT WIDB_SiteSub_dStem_TT
clear matchingrows WIDB_SiteSub_dStem_TT_temp WIDB_SiteSub_dSoil_TT_temp

%Match site name
matching = ismember(WIDB_Table.Site_Name, siteid.site_long_name{u});
WIDB_SiteSub_Table = WIDB_Table(matching, :);

%Separate into soil-water and plant-water samples
WIDB_SiteSub_Soil_Table = WIDB_SiteSub_Table(strcmp(WIDB_SiteSub_Table.Type, 'Soil'), :);
WIDB_SiteSub_Stem_Table = WIDB_SiteSub_Table(strcmp(WIDB_SiteSub_Table.Type, 'Stem'), :);

WIDB_SiteSub_dSoil_TT   = table2timetable(WIDB_SiteSub_Soil_Table);
WIDB_SiteSub_dStem_TT   = table2timetable(WIDB_SiteSub_Stem_Table);

%Match sample date
matchingrows               = find(WIDB_SiteSub_dStem_TT.Collection_Date == dateshift(SampDate.Start_Sample_Day(u),'start','day'));
WIDB_SiteSub_dStem_TT_temp = WIDB_SiteSub_dStem_TT(matchingrows,:);

%Calculate and save the mean, number of samples and std. dev
dT_18O_summ_out(u,1) = mean(WIDB_SiteSub_dStem_TT_temp.d18O,'all','omitnan');
dT_18O_summ_out(u,3) = nnz(~isnan(WIDB_SiteSub_dStem_TT_temp.d18O));
dT_18O_summ_out(u,2) = std(WIDB_SiteSub_dStem_TT_temp.d18O,0,1,"omitnan")/sqrt(dT_18O_summ_out(u,3));

dT_2H_summ_out(u,1) = mean(WIDB_SiteSub_dStem_TT_temp.d2H,'all','omitnan');
dT_2H_summ_out(u,3) =nnz(~isnan(WIDB_SiteSub_dStem_TT_temp.d2H));
dT_2H_summ_out(u,2) =std(WIDB_SiteSub_dStem_TT_temp.d2H,0,1,"omitnan")/sqrt(dT_2H_summ_out(u,3));

cd(fullfile(code_path))
writetimetable(WIDB_SiteSub_dStem_TT_temp,'dT_Matlab.xlsx','Sheet',siteid.site_id{u});

%dS
clear matchingrows
matchingrows = find(WIDB_SiteSub_dSoil_TT.Collection_Date == dateshift(SampDate.Start_Sample_Day(u),'start','day') & WIDB_SiteSub_dSoil_TT.Depth_meters <= 0.25);
WIDB_SiteSub_dSoil_TT_temp      =     WIDB_SiteSub_dSoil_TT(matchingrows,:);

%Calculate and save the mean, number of samples and std. dev
dSoil_18O_summ_out(u,1) = mean(WIDB_SiteSub_dSoil_TT_temp.d18O,'all','omitnan');
dSoil_18O_summ_out(u,3) = nnz(~isnan(WIDB_SiteSub_dSoil_TT_temp.d18O));
dSoil_18O_summ_out(u,2) = std(WIDB_SiteSub_dSoil_TT_temp.d18O,0,1,"omitnan")/sqrt(dSoil_18O_summ_out(u,3));

dSoil_2H_summ_out(u,1) = mean(WIDB_SiteSub_dSoil_TT_temp.d2H,'all','omitnan');
dSoil_2H_summ_out(u,3) =nnz(~isnan(WIDB_SiteSub_dSoil_TT_temp.d2H));
dSoil_2H_summ_out(u,2) =std(WIDB_SiteSub_dSoil_TT_temp.d2H,0,1,"omitnan")/sqrt(dSoil_2H_summ_out(u,3));

cd(fullfile(code_path))
writetimetable(WIDB_SiteSub_dSoil_TT_temp,'dS_Matlab.xlsx','Sheet',siteid.site_id{u});

end

%Save outputs
cd(fullfile(code_path))

writematrix(dT_18O_summ_out,'dT_Matlab.xlsx','Sheet','dT 18O Summ');
writematrix(dT_2H_summ_out,'dT_Matlab.xlsx','Sheet','dT 2H Summ');
writematrix(dT_18O_summ_out,'dT_18O_summ_out.xlsx','Sheet','dT 18O Summ');
writematrix(dT_2H_summ_out,'dT_2H_summ_out.xlsx','Sheet','dT 2H Summ');

writematrix(dSoil_18O_summ_out,'dS_Matlab.xlsx','Sheet','dS 18O Summ');
writematrix(dSoil_2H_summ_out,'dS_Matlab.xlsx','Sheet','dS 2H Summ');
writematrix(dSoil_18O_summ_out,'dS_18O_summ_out.xlsx','Sheet','dS 18O Summ');
writematrix(dSoil_2H_summ_out,'dS_2H_summ_out.xlsx','Sheet','dS 2H Summ');


%% δE(α*, ε*, δS, hA, n, εk) and T/ET: 3-day 
%Description: Isotopic composition of evaporated soil-water, Monte Carlo Analysis and Transpiration Fraction (fT; T/ET)

%Load Variables
cd (code_path)
Ts_summ_out = readmatrix('Ts_summ_out_Matlab.xlsx');
Ta_summ_out = readmatrix('Ta_summ_out_Matlab.xlsx');
RH_summ_out = readmatrix('RH_summ_out_Matlab.xlsx');
dSoil_18O_summ_out = readmatrix('dS_18O_summ_out.xlsx');
dSoil_2H_summ_out  = readmatrix('dS_2H_summ_out.xlsx');
dA_18O_summ_out    = readmatrix('dA_18O_Lvl_summ_out.xlsx');
dA_2H_summ_out     = readmatrix('dA_2H_Lvl_summ_out.xlsx');

%% ET

for sampdate =  1:56

%Initiate empty arrays    
Ts = (zeros(1,nn))'; Ta = (zeros(1,nn))'; RH = (zeros(1,nn))'; hA = (zeros(1,nn))'; n = (zeros(1,nn))';

dS_18O = (zeros(1,nn))'; dS_2H = (zeros(1,nn))'; dA_18O = (zeros(1,nn))'; dA_2H = (zeros(1,nn))';
dE_18O = (zeros(1,nn))'; dE_2H = (zeros(1,nn))';


for i = 1:1000

%Sample variable from distribution with mean and std. dev of calculated above   
Ts_temp       = normrnd(Ts_summ_out(sampdate,1),Ts_summ_out(sampdate,2)*sqrt(Ts_summ_out(sampdate,3))); Ts(i) = Ts_temp;
Ta_temp       = normrnd(Ta_summ_out(sampdate,1),Ta_summ_out(sampdate,2)*sqrt(Ta_summ_out(sampdate,3))); Ta(i) = Ta_temp;
RH_temp       = normrnd(RH_summ_out(sampdate,1),RH_summ_out(sampdate,2)*sqrt(RH_summ_out(sampdate,3))); RH(i) = RH_temp;

dS_18O_temp   = normrnd(dSoil_18O_summ_out(sampdate,1),dSoil_18O_summ_out(sampdate,2)*sqrt(dSoil_18O_summ_out(sampdate,3)));dS_18O (i) = dS_18O_temp;
dS_2H_temp    = normrnd(dSoil_2H_summ_out(sampdate,1),dSoil_2H_summ_out(sampdate,2)*sqrt(dSoil_2H_summ_out(sampdate,3)));dS_2H (i) = dS_2H_temp;
dA_18O_temp   = normrnd(dA_18O_summ_out(sampdate,1),dA_18O_summ_out(sampdate,2)*sqrt(dA_18O_summ_out(sampdate,3)));dA_18O (i) = dA_18O_temp;
dA_2H_temp    = normrnd(dA_2H_summ_out(sampdate,1),dA_2H_summ_out(sampdate,2)*sqrt(dA_2H_summ_out(sampdate,3)));dA_2H (i) = dA_2H_temp;

% Variable (2): α*
alpha_18O_temp = 1/(exp((-7.685 + 6.7123.*(10^3./(Ts_temp+273.15)) - 1.6664.*(10^6./(Ts_temp+273.15).^2) + 0.35041.*(10^9./(Ts_temp+273.15).^3))/10^3));
alpha_2H_temp  = 1/(exp((1158.8.*((Ts_temp+273.15).^3./10^9) - 1620.1.*((Ts_temp+273.15).^2./10^6) + 794.84.*((Ts_temp+273.15)./10^3)- 161.04 + 2.9992.*(10^9./(Ts_temp+273.15).^3))/10^3));

% Variable (3): ε*
eps_eq_18O_temp = (1-alpha_18O_temp)*10^3;
eps_eq_2H_temp  = (1-alpha_2H_temp)*10^3;

% Variable (5): hA 
esat_air  = (6.11*10.^((7.5.*Ta_temp)./(237.3+Ta_temp)))./10;
esat_soil = (6.11*10.^((7.5.*Ts_temp)./(237.3+Ts_temp)))./10;
hA_temp= (RH_temp*esat_air/esat_soil)/100; hA (i) = hA_temp;

if hA_temp > 0.8
   hA_temp = 0.8;
end

hA (nn) = hA_temp;

% Variable (6): n 
%Sample from normal distribution over possible range of "n"
n_temp = 0.5 + (1-0.5)*rand(1,1); n(i) = n_temp;

% Variable (7): εk
ek_18O_temp = n_temp*(1 - hA_temp)*(1-D_18O)*1000;
ek_2H_temp  = n_temp*(1 - hA_temp)*(1-D_2H)*1000;

%Craig and Gordon (1965)
dE_18O_temp = (((alpha_18O_temp*dS_18O_temp-hA_temp*dA_18O_temp))-(eps_eq_18O_temp+ek_18O_temp)/((1-hA_temp)+(ek_18O_temp/1000))); dE_18O (i) = dE_18O_temp;
dE_2H_temp  = (((alpha_2H_temp*dS_2H_temp-hA_temp*dA_2H_temp))-(eps_eq_2H_temp+ek_18O_temp)/((1-hA_temp)+(ek_2H_temp/1000)));  dE_2H (i) = dE_2H_temp;

end

output_concat = horzcat(Ts, Ta, RH, hA, n, dS_18O, dS_2H, dA_18O, dA_2H, dE_18O, dE_2H);
output_concat_table = array2table(output_concat,"VariableNames",{'Ts', 'Ta', 'RH', 'hA', 'n', 'dS_18O', 'dS_2H', 'dA_18O', 'dA_2H','dE_18O', 'dE_2H'});
writetable(output_concat_table,'All-Vars_MC-Output_Matlab.xlsx','Sheet',siteid.site_id{sampdate});

output_sum_mean(sampdate,:)    = mean(output_concat);
output_sum_std(sampdate,:)     = std(output_concat);

writematrix(output_sum_mean,'summary_dE.xlsx','Sheet','mean');
writematrix(output_sum_std,'summary_dE.xlsx','Sheet','std');

end

%% Pull in dET values for time windows & calculate composite (mean) dET and composite (total) standard error

%End members together 
dET_18O_1day_summ_out = readtable('Keeling_Summary_Outputs_tw1day.xlsx','Sheet','18O');
dET_2H_1day_summ_out  = readtable('Keeling_Summary_Outputs_tw1day.xlsx','Sheet','2H');

dET_18O_3day_summ_out = readtable('Keeling_Summary_Outputs_tw3day.xlsx','Sheet','18O');
dET_2H_3day_summ_out  = readtable('Keeling_Summary_Outputs_tw3day.xlsx','Sheet','2H');

dET_18O_5day_summ_out = readtable('Keeling_Summary_Outputs_tw5day.xlsx','Sheet','18O');
dET_2H_5day_summ_out  = readtable('Keeling_Summary_Outputs_tw5day.xlsx','Sheet','2H');

dET_18O_7day_summ_out = readtable('Keeling_Summary_Outputs_tw7day.xlsx','Sheet','18O');
dET_2H_7day_summ_out  = readtable('Keeling_Summary_Outputs_tw7day.xlsx','Sheet','2H');

dETval_18O_tw     = table(dET_18O_1day_summ_out.LinSlope,dET_18O_3day_summ_out.LinSlope,dET_18O_5day_summ_out.LinSlope,dET_18O_7day_summ_out.LinSlope,'VariableNames',{'dET18O_1day','dET18O_3day','dET18O_5day','dET18O_7day'});
dETval_2H_tw      = table(dET_2H_1day_summ_out.LinSlope,dET_2H_3day_summ_out.LinSlope,dET_2H_5day_summ_out.LinSlope,dET_2H_7day_summ_out.LinSlope,'VariableNames',{'dET2H_1day','dET2H_3day','dET2H_5day','dET2H_7day'});

dETSE_18O_tw      = table(dET_18O_1day_summ_out.LinSE,dET_18O_3day_summ_out.LinSE,dET_18O_5day_summ_out.LinSE,dET_18O_7day_summ_out.LinSE,'VariableNames',{'dET18O_1day','dET18O_3day','dET18O_5day','dET18O_7day'});
dETSE_2H_tw       = table(dET_2H_1day_summ_out.LinSE,dET_2H_3day_summ_out.LinSE,dET_2H_5day_summ_out.LinSE,dET_2H_7day_summ_out.LinSE,'VariableNames',{'dET2H_1day','dET2H_3day','dET2H_5day','dET2H_7day'});

dETval_comp_18O   = mean(dETval_18O_tw,2,"omitnan");
dETval_comp_2H    = mean(dETval_2H_tw,2,"omitnan");

dET_SEpooled_18O   = sqrt((dETSE_18O_tw.dET18O_1day.^2)+(dETSE_18O_tw.dET18O_3day.^2)+(dETSE_18O_tw.dET18O_5day.^2)+(dETSE_18O_tw.dET18O_7day.^2))/sqrt(size(dETSE_18O_tw,2));
dET_SEvarMeans_18O = std(dETval_18O_tw,0,2)./sqrt(size(dETSE_18O_tw,2));
dET_SEtotal_18O    = sqrt(dET_SEpooled_18O.^2+dET_SEvarMeans_18O.^2);

%ONAQ3 only has 3 dET values, not dET 1-day
dET_SEpooled_18O(49,1)   = sqrt((dETSE_18O_tw.dET18O_3day(49).^2)+(dETSE_18O_tw.dET18O_5day(49).^2)+(dETSE_18O_tw.dET18O_7day(49).^2))/sqrt(3);
dET_SEvarMeans_18O(49,1) = std(dETval_18O_tw(49,2:4),0,2)./sqrt(3);
dET_SEtotal_18O.std(49)    = sqrt(dET_SEpooled_18O(49)^2+dET_SEvarMeans_18O.std(49)^2);

dET_SEpooled_2H    = sqrt((dETSE_2H_tw.dET2H_1day.^2)+(dETSE_2H_tw.dET2H_3day.^2)+(dETSE_2H_tw.dET2H_5day.^2)+(dETSE_2H_tw.dET2H_7day.^2))/sqrt(size(dETSE_2H_tw,2));
dET_SEvarMeans_2H  = std(dETval_2H_tw,0,2)./sqrt(size(dETSE_2H_tw,2));
dET_SEtotal_2H     = sqrt(dET_SEpooled_2H.^2+dET_SEvarMeans_2H.^2);

%ONAQ3 only has 3 dET values, not dET 1-day
dET_SEpooled_2H(49,1)   = sqrt((dETSE_2H_tw.dET2H_3day(49).^2)+(dETSE_2H_tw.dET2H_5day(49).^2)+(dETSE_2H_tw.dET2H_7day(49).^2))/sqrt(3);
dET_SEvarMeans_2H(49,1) = std(dETval_2H_tw(49,2:4),0,2)./sqrt(3);
dET_SEtotal_2H.std(49)    = sqrt(dET_SEpooled_2H(49)^2+dET_SEvarMeans_2H.std(49)^2);

% Composite (mean) dET 18O --> dETval_comp_18O 
% Composite (mean) dET 2H  --> dETval_comp_2H

% Composite (total) SE dET 18O --> dET_SEtotal_18O 
% Composite (total) SE dET 2H  --> dET_SEtotal_2H


%% Pull in dE values and standard error from Monte Carlo 

dE_mean = readtable('summary_dE.xlsx','Sheet','mean');
dE_SE   = readtable('summary_dE.xlsx','Sheet','std');

dEval_18O = dE_mean.Var10;
dEval_2H  = dE_mean.Var11;

dE_SE_18O = dE_SE.Var10;
dE_SE_2H  = dE_SE.Var11;

% dE 18O mean --> dEval_18O
% dE 2H  mean --> dEval_2H

% dE 18O SE --> dE_SE_18O
% dE 2H SE  --> dE_SE_2H

%% Pull in dT values and standard error from spectral corrections

dT_18O    = readmatrix('dT_18O_summ_out_CRDSUpdate.xlsx');
dT_2H     = readmatrix('dT_2H_summ_out_CRDSUpdate.xlsx');

dTval_18O = dT_18O(:,1);
dTval_2H  = dT_2H(:,1);

dT_SE_18O = dT_18O(:,2);
dT_SE_2H  = dT_2H(:,2);


% dT 18O mean --> dTval_18O
% dT 2H  mean --> dTval_2H

% dT 18O SE --> dT_SE_18O
% dT 2H SE  --> dT_SE_2H

%% Calculate T/ET at the Sample Date Level

% T/ET and SE WITHOUT BIAS
fTval_18O = (dETval_comp_18O.mean - dEval_18O)./(dTval_18O - dEval_18O);
fTval_2H = (dETval_comp_2H.mean - dEval_2H)./(dTval_2H - dEval_2H);

%Gaussian Error Propagation (Pre-Bias)
%Phillips & Gregg, 2001 Eq. A4

%First Part  = ((1/(dT-dE))^2)*(SE(dET))^2
%Second Part = ((fT/(dT-dE))^2)*(SE(dT))^2
%Third Part  = (((fT-1)/(dT-dE))^2)*(SE(dE))^2

%18O
first_part_SE_18O  = ((1./(dTval_18O-dEval_18O)).^2).*(dET_SEtotal_18O.std).^2;
second_part_SE_18O = ((-fTval_18O./(dTval_18O-dEval_18O)).^2).*(dT_SE_18O).^2;
third_part_SE_18O  = (((fTval_18O - 1)./(dTval_18O-dEval_18O)).^2).*(dE_SE_18O).^2;

fT_SE_18O = sqrt(first_part_SE_18O+second_part_SE_18O+third_part_SE_18O);

%2H
first_part_SE_2H  = ((1./(dTval_2H-dEval_2H)).^2).*(dET_SEtotal_2H.std).^2;
second_part_SE_2H = ((-fTval_2H./(dTval_2H-dEval_2H)).^2).*(dT_SE_2H).^2;
third_part_SE_2H  = (((fTval_2H - 1)./(dTval_2H-dEval_2H)).^2).*(dE_SE_2H).^2;

fT_SE_2H = sqrt(first_part_SE_2H+second_part_SE_2H+third_part_SE_2H);

%% T/ET CVE BIAS 
%BIAS Equation: Allen & Kirchner, 2022 Eq. 4
%BIAS Values from Wen et al., 2023: δ18O: -0.79±0.34
%BIAS Values from Wen et al., 2023:  δ2H: -8.21±0.67

%Three part equation but BIAS(dET) = 0 & BIAS(dE) = 0
%BIAS (fT) = -fT*BIAS(dT)/(dT-dE)

%18O 
BIAS_dTval_18O = 0.79;
BIAS_dT_SE_18O = 0.34;

BIAS_fT_18O      = (-fTval_18O.*BIAS_dTval_18O)./(dTval_18O-dEval_18O);
fTval_18O_biased = fTval_18O + BIAS_fT_18O;

%2H
BIAS_dTval_2H = 8.21;
BIAS_dT_SE_2H = 0.67;

BIAS_fT_2H      = (-fTval_2H.*BIAS_dTval_2H)./(dTval_2H-dEval_2H);
fTval_2H_biased = fTval_2H + BIAS_fT_2H;

%% T/ET CVE BIAS SE

%BIAS Equation: Allen & Kirchner, 2022 Eq. 10
%BIAS Values from Wen et al., 2023: δ18O: -0.79±0.34
%BIAS Values from Wen et al., 2023:  δ2H: -8.21±0.67

%First Part = 1/abs(δT-δE)
%Second Part = var(δET)
%Third Part = (fT^2)*var(dT)
%Fourth Part = ((1-fT)^2)*(var(δE))
%Fifth Part = 0; BIAS(dET) = 0;
%Sixth Part = (fT^2)*SE_BIAS(δT)^2
%Seventh Part = 0; BIAS(dE) = 0;

%18O 
first_part_biased_SE_18O  = 1./abs(dTval_18O-dEval_18O);
second_part_biased_SE_18O = (dET_SEtotal_18O).^2;
third_part_biased_SE_18O  = (fTval_18O_biased.^2).*(dT_SE_18O.^2);
fourth_part_biased_SE_18O = ((1-fTval_18O_biased).^2).*(dE_SE_18O.^2);
sixth_part_biased_SE_18O  = (fTval_18O_biased.^2).*(BIAS_dT_SE_18O.^2);
fT_SE_18O_biased = first_part_biased_SE_18O.*sqrt(first_part_biased_SE_18O+second_part_biased_SE_18O+third_part_biased_SE_18O+fourth_part_biased_SE_18O+sixth_part_biased_SE_18O+sixth_part_biased_SE_18O);

%2H 
first_part_biased_SE_2H  = 1./abs(dTval_2H-dEval_2H);
second_part_biased_SE_2H = (dET_SEtotal_2H).^2;
third_part_biased_SE_2H  = (fTval_2H_biased.^2).*(dT_SE_2H.^2);
fourth_part_biased_SE_2H = ((1-fTval_2H_biased).^2).*(dE_SE_2H.^2);
sixth_part_biased_SE_2H  = (fTval_2H_biased.^2).*(BIAS_dT_SE_2H.^2);
fT_SE_2H_biased = first_part_biased_SE_2H.*sqrt(first_part_biased_SE_2H+second_part_biased_SE_2H+third_part_biased_SE_2H+fourth_part_biased_SE_2H+sixth_part_biased_SE_2H+sixth_part_biased_SE_2H);

%% TABLES 

% TABLE 2: By Sample Date
Table2 = table(siteid.site_id,fTval_18O_biased,fT_SE_18O_biased.std,fTval_2H_biased,fT_SE_2H_biased.std,'VariableNames',{'SiteID','fT_biased_18O','fT_biased_SE_18O','fT_biased_2H','fT_biased_SE_2H'});

% Table 3: By Site
Table2a = table(siteid.site_name,fTval_18O_biased,fT_SE_18O_biased.std,fTval_2H_biased,fT_SE_2H_biased.std,'VariableNames',{'SiteName','fT_biased_18O','fT_biased_SE_18O','fT_biased_2H','fT_biased_SE_2H'});

Table2a.SiteName = categorical(Table2a.SiteName);

Site.empty = table(1);

[C,ia] = unique(siteid.site_name);

%Grab each site
for i = 1:length(ia)

Site.(siteid.site_name{ia(i)}) =Table2a(Table2a.SiteName==siteid.site_name{ia(i)},:);

end


%Calculations for each site
for i = 1:length(ia)

%Mean fT by Site
site_mean_fT_18O(i,1) = mean(Site.(siteid.site_name{ia(i)}).fT_biased_18O);
site_mean_fT_2H(i,1)  = mean(Site.(siteid.site_name{ia(i)}).fT_biased_2H);

%SE by site
SE_pooled_18O (i,1) = sqrt(sum((Site.(siteid.site_name{ia(i)}).fT_biased_SE_18O).^2)/numel(Site.(siteid.site_name{ia(i)}).fT_biased_SE_18O));
SE_of_mean_18O(i,1) = std(Site.(siteid.site_name{ia(i)}).fT_biased_18O)/sqrt(numel(Site.(siteid.site_name{ia(i)}).fT_biased_SE_18O));
site_SE_18O (i,1)   = sqrt(SE_pooled_18O(i,1)^2+SE_of_mean_18O(i,1)^2);

SE_pooled_2H (i,1) = sqrt(sum((Site.(siteid.site_name{ia(i)}).fT_biased_SE_2H).^2)/numel(Site.(siteid.site_name{ia(i)}).fT_biased_SE_2H));
SE_of_mean_2H(i,1) = std(Site.(siteid.site_name{ia(i)}).fT_biased_2H)/sqrt(numel(Site.(siteid.site_name{ia(i)}).fT_biased_SE_2H));
site_SE_2H (i,1)   = sqrt(SE_pooled_2H(i,1)^2+SE_of_mean_2H(i,1)^2);

end

Table3 = table(C,site_mean_fT_18O,site_SE_18O,site_mean_fT_2H,site_SE_2H,'VariableNames',{'Site','fT_biased_18O','fT_biased_SE_18O','fT_biased_2H','fT_biased_SE_2H'});


% Table 4: By Ecosystem Type
Table2b = table(siteid.ecosystemType,fTval_18O_biased,fT_SE_18O_biased.std,fTval_2H_biased,fT_SE_2H_biased.std,'VariableNames',{'EcosystemType','fT_biased_18O','fT_biased_SE_18O','fT_biased_2H','fT_biased_SE_2H'});

Table2b.EcosystemType = categorical(Table2b.EcosystemType);

Ecosystem.empty = table(1);

[D,ib] = unique(siteid.ecosystemType);

%Grab each EcosystemType
for i = 1:length(ib)

Ecosystem.(siteid.ecosystemType{ib(i)}) =Table2b(Table2b.EcosystemType==siteid.ecosystemType{ib(i)},:);

end


%Calculations for each ecosystem type
for i = 1:length(ib)

%Mean fT by EcoType
eco_mean_fT_18O(i,1) = mean(Ecosystem.(siteid.ecosystemType{ib(i)}).fT_biased_18O);
eco_mean_fT_2H(i,1)  = mean(Ecosystem.(siteid.ecosystemType{ib(i)}).fT_biased_2H);

%SE by EcoType
SE_pooled_18O_eco (i,1) = sqrt(sum((Ecosystem.(siteid.ecosystemType{ib(i)}).fT_biased_SE_18O).^2)/numel(Ecosystem.(siteid.ecosystemType{ib(i)}).fT_biased_SE_18O));
SE_of_mean_18O_eco(i,1) = std(Ecosystem.(siteid.ecosystemType{ib(i)}).fT_biased_18O)/sqrt(numel(Ecosystem.(siteid.ecosystemType{ib(i)}).fT_biased_SE_18O));
eco_SE_18O (i,1)    = sqrt(SE_pooled_18O_eco(i,1)^2+SE_of_mean_18O_eco(i,1)^2);

SE_pooled_2H_eco (i,1) = sqrt(sum((Ecosystem.(siteid.ecosystemType{ib(i)}).fT_biased_SE_2H).^2)/numel(Ecosystem.(siteid.ecosystemType{ib(i)}).fT_biased_SE_2H));
SE_of_mean_2H_eco(i,1) = std(Ecosystem.(siteid.ecosystemType{ib(i)}).fT_biased_2H)/sqrt(numel(Ecosystem.(siteid.ecosystemType{ib(i)}).fT_biased_SE_2H));
eco_SE_2H (i,1)    = sqrt(SE_pooled_2H_eco(i,1)^2+SE_of_mean_2H_eco(i,1)^2);

end

Table4 = table(D,eco_mean_fT_18O,eco_SE_18O,eco_mean_fT_2H,eco_SE_2H,'VariableNames',{'EcosystemType','fT_biased_18O','fT_biased_SE_18O','fT_biased_2H','fT_biased_SE_2H'});


% Table 5: By Season Type
Table2c = table(siteid.season,fTval_18O_biased,fT_SE_18O_biased.std,fTval_2H_biased,fT_SE_2H_biased.std,'VariableNames',{'Season','fT_biased_18O','fT_biased_SE_18O','fT_biased_2H','fT_biased_SE_2H'});

Table2c.Season = categorical(Table2c.Season);

Season.empty = table(1);

[E,ic] = unique(siteid.season);

%Grab each Season
Season.Fall   =Table2c(Table2c.Season=="Fall",:);
Season.Spring =Table2c(Table2c.Season=="Spring",:);
Season.Summer =Table2c(Table2c.Season=="Summer",:);

%Calculations for each Season
for i = 1:3

%Mean fT by Season
season_mean_fT_18O(i,1) = mean(Season.(siteid.season{ic(i)}).fT_biased_18O);
season_mean_fT_2H(i,1) = mean(Season.(siteid.season{ic(i)}).fT_biased_2H);

%SE by Season
SE_pooled_18O_season (i,1) = sqrt(sum((Season.(siteid.season{ic(i)}).fT_biased_SE_18O).^2)/numel(Season.(siteid.season{ic(i)}).fT_biased_SE_18O));
SE_of_mean_18O_season(i,1) = std(Season.(siteid.season{ic(i)}).fT_biased_18O)/sqrt(numel(Season.(siteid.season{ic(i)}).fT_biased_SE_18O));
season_SE_18O (i,1)        = sqrt(SE_pooled_18O_season(i,1)^2+SE_of_mean_18O_season(i,1)^2);

SE_pooled_2H_season (i,1) = sqrt(sum((Season.(siteid.season{ic(i)}).fT_biased_SE_2H).^2)/numel(Season.(siteid.season{ic(i)}).fT_biased_SE_2H));
SE_of_mean_2H_season(i,1) = std(Season.(siteid.season{ic(i)}).fT_biased_2H)/sqrt(numel(Season.(siteid.season{ic(i)}).fT_biased_SE_2H));
season_SE_2H (i,1)        = sqrt(SE_pooled_2H_season(i,1)^2+SE_of_mean_2H_season(i,1)^2);

end

Table5 = table(E,season_mean_fT_18O,season_SE_18O,season_mean_fT_2H,season_SE_2H,'VariableNames',{'Season','fT_biased_18O','fT_biased_SE_18O','fT_biased_2H','fT_biased_SE_2H'});



%% Pull in LAI, NDVI, and ET

NEONSitesMCD15A2H061results  = readtable('NEON_Appears_Raw_Data.xlsx','Sheet','NEON_Appears_Raw_Data_05132024');
NEONSitesMOD13A1061results   = readtable('NEON_Appears_Raw_Data.xlsx','Sheet','NEONSitesMOD13A1061results');
NEONSitesMOD16A2GF061results = readtable('NEON_Appears_Raw_Data.xlsx','Sheet','NEONSitesMOD16A2GF061results');


%%

Table2d = table(siteid.site_name, siteid.site_id,siteid.ecosystemType,siteid.season,SampDate.Start_Sample_Day,fTval_18O_biased,fT_SE_18O_biased.std,fTval_2H_biased,fT_SE_2H_biased.std,'VariableNames',{'SiteName','SiteID','EcosystemType','Season','SampleDate','fT_biased_18O','fT_biased_SE_18O','fT_biased_2H','fT_biased_SE_2H'});
%
store = table(siteid.site_name,'VariableNames',{'SiteName'});

store.LAI  = NaN(56,1);
store.NDVI = NaN(56,1);
store.ET_GapFilled = NaN(56,1);
store.ET_EC = NaN(56,1);

for u = [1;4;10;14;16;19;25;27;31;36;41;44;47]'

%LAI
clear matching LAI_TableSubset LAI_TT
matching = ismember(NEONSitesMCD15A2H061results.ID,siteid.site_name{u});
LAI_TableSubset = NEONSitesMCD15A2H061results(matching, :);
NEON.(siteid.site_name{u}).LAI = table2timetable(LAI_TableSubset);
NEON.(siteid.site_name{u}).LAI = sortrows(NEON.(siteid.site_name{u}).LAI);
NEON.(siteid.site_name{u}).LAI = unique(NEON.(siteid.site_name{u}).LAI);

%LAI Interpolate
LAI_TT = timetable(NEON.(siteid.site_name{u}).LAI.Date,(NEON.(siteid.site_name{u}).LAI.MCD15A2H_061_Lai_500m),'VariableNames',{'LAI'});
LAI_TT = sortrows(LAI_TT);
LAI_TT = unique(LAI_TT);
NEON.(siteid.site_name{u}).LAI_TT_daily = retime(LAI_TT,'daily','spline');


%LAI values on Sample Dates
clear matching LAI_interp_subset
matching          = ismember(dateshift(NEON.(siteid.site_name{u}).LAI_TT_daily.Time,'start','day'), dateshift(Table2d.SampleDate(ismember(Table2d.SiteName, siteid.site_name{u})),'start','day'));
LAI_interp_subset = NEON.(siteid.site_name{u}).LAI_TT_daily(matching, :);
store.LAI(ismember(store.SiteName, siteid.site_name{u})) = LAI_interp_subset.LAI;

%NDVI
clear matching EVI_TableSubset
matching = ismember(NEONSitesMOD13A1061results.ID,siteid.site_name{u});
EVI_TableSubset = NEONSitesMOD13A1061results(matching, :);
NEON.(siteid.site_name{u}).EVI = table2timetable(EVI_TableSubset);
NEON.(siteid.site_name{u}).EVI = sortrows(NEON.(siteid.site_name{u}).EVI);
NEON.(siteid.site_name{u}).EVI = unique(NEON.(siteid.site_name{u}).EVI);

%NDVI Interpolate
clear EVI_TT
EVI_TT = timetable(NEON.(siteid.site_name{u}).EVI.Date,(NEON.(siteid.site_name{u}).EVI.MOD13A1_061__500m_16_days_NDVI),'VariableNames',{'NDVI'});
EVI_TT = sortrows(EVI_TT);
EVI_TT = unique(EVI_TT);
NEON.(siteid.site_name{u}).EVI_TT_daily = retime(EVI_TT,'daily','spline');

%NDVI values on Sample Dates
clear matching EVI_interp_subset
matching          = ismember(dateshift(NEON.(siteid.site_name{u}).EVI_TT_daily.Time,'start','day'), dateshift(Table2d.SampleDate(ismember(Table2d.SiteName, siteid.site_name{u})),'start','day'));
EVI_interp_subset = NEON.(siteid.site_name{u}).EVI_TT_daily(matching, :);
store.NDVI(ismember(store.SiteName, siteid.site_name{u})) = EVI_interp_subset.NDVI;

%ET
clear matching
matching = ismember(NEONSitesMOD16A2GF061results.ID,siteid.site_name{u});
ET_TableSubset = NEONSitesMOD16A2GF061results(matching, :);
NEON.(siteid.site_name{u}).ET = table2timetable(ET_TableSubset);
NEON.(siteid.site_name{u}).ET = sortrows(NEON.(siteid.site_name{u}).ET);
NEON.(siteid.site_name{u}).ET = unique(NEON.(siteid.site_name{u}).ET);

%NIWO and WREF don't have Gap Filled ET
if u == 14 || u == 41

    continue
end

%ET Gap Filled values on sample dates
clear matching ET_interp_subset
matching = ismember(dateshift(ONEFLUX.(siteid.site_name{u}).Flux15_DD_TT.Time,'start','day'), dateshift(Table2d.SampleDate(ismember(Table2d.SiteName, siteid.site_name{u})),'start','day'));
ET_interp_subset = ONEFLUX.(siteid.site_name{u}).Flux15_DD_TT.ET(matching, :);

store.ET_GapFilled(ismember(store.SiteName, siteid.site_name{u})) = ET_interp_subset;

%ET EC
clear ET_EC matching ET_EC_sub
ET_EC     = readtable('ET_EC_Sites.xlsx','Sheet',siteid.site_name{u});

if u == 19 || u == 27 || u == 47

ET_EC.ET_nsae_mmday(cellfun(@isempty,ET_EC.ET_nsae_mmday)) = {NaN};
matching  = ismember(dateshift(ET_EC.Datetime,'start','day'), dateshift(Table2d.SampleDate(ismember(Table2d.SiteName, siteid.site_name{u})),'start','day'));
ET_EC_sub = ET_EC.ET_nsae_mmday(matching);
store.ET_EC(ismember(store.SiteName, siteid.site_name{u})) = str2double(ET_EC_sub);

else


matching  = ismember(dateshift(ET_EC.Datetime,'start','day'), dateshift(Table2d.SampleDate(ismember(Table2d.SiteName, siteid.site_name{u})),'start','day'));
ET_EC_sub = ET_EC.ET_nsae_mmday(matching);
store.ET_EC(ismember(store.SiteName, siteid.site_name{u})) = ET_EC_sub;

if u == 41
clear ET_EC matching ET_EC_sub
ET_EC     = readtable('ET_EC_Sites.xlsx','Sheet',siteid.site_name{u});
matching  = ismember(dateshift(ET_EC.Datetime,'start','day'), dateshift(Table2d.SampleDate(ismember(Table2d.SiteName, siteid.site_name{u})),'start','day'));
ET_EC_sub = ET_EC.ET_turb_mmday(matching);
store.ET_EC(ismember(store.SiteName, siteid.site_name{u})) = ET_EC_sub;

end
end
end

%CORRECTION: WREF3 LAI = WREF2 LAI 
store.LAI(43) = store.LAI(42);


store.ET_Combined = store.ET_EC;

for u = [4; 11; 16; 19; 20; 21; 28; 29; 30; 39; 44; 45; 49; 50; 51; 56]'

store.ET_Combined(u) = store.ET_GapFilled(u);

end


Table2D = [Table2d,store(:,2:end)];
Table2D = removevars(Table2D,["ET_GapFilled","ET_EC"]);

TranFlux_18O = NaN(56,1);
TranFlux_2H  = NaN(56,1);

EFlux_18O = NaN(56,1);
EFlux_2H  = NaN(56,1);

for i = 1:56

if Table2D.fT_biased_18O(i)>1
   TranFlux_18O(i) = 1*Table2D.ET_Combined(i);
   EFlux_18O (i)   = 0;
elseif Table2D.fT_biased_18O(i)<0
    TranFlux_18O(i) = 0*Table2D.ET_Combined(i);
    EFlux_18O (i)   = 1*Table2D.ET_Combined(i);
else
    TranFlux_18O(i) = Table2D.fT_biased_18O(i)*Table2D.ET_Combined(i);
    EFlux_18O(i) = (1-Table2D.fT_biased_18O(i))*Table2D.ET_Combined(i);
end

if Table2D.fT_biased_2H(i)>1
   TranFlux_2H(i) = 1*Table2D.ET_Combined(i);
   EFlux_2H (i)   = 0;
elseif Table2D.fT_biased_2H(i)<0
    TranFlux_2H(i) = 0*Table2D.ET_Combined(i);
    EFlux_2H (i)   = 1*Table2D.ET_Combined(i);
else
    TranFlux_2H(i) = Table2D.fT_biased_2H(i)*Table2D.ET_Combined(i);
    EFlux_2H(i) = (1-Table2D.fT_biased_2H(i))*Table2D.ET_Combined(i);
end

end

Table7 = table(siteid.site_id,SampDate.Start_Sample_Day,store.ET_EC,store.ET_GapFilled,store.ET_Combined,Table2D.fT_biased_18O,TranFlux_18O,Table2D.fT_biased_2H,TranFlux_2H,Table2D.NDVI,Table2D.LAI,'VariableNames',{'SiteID','SampleDate','EC_EC','ET_GapFilled','ET_Combined','fT_biased_18O','TFlux_18O','fT_biased_2H','TFlux_2H','NDVI','LAI'});
EFlux = table(EFlux_18O,EFlux_2H,'VariableNames',{'EFlux18O','EFlux2H'});

%Correlations Table6
%Table6: All Sample Dates
Corr = nan(12,3);

%T/ET (δ18O) vs LAI
row = 1;
[Corr(row,1),Corr(row,2)] = corr(Table7.LAI, Table7.fT_biased_18O,"Type","Kendall");
lm = fitlm(Table7.LAI, Table7.fT_biased_18O);
Corr(row,3) = lm.Rsquared.Ordinary;

%T/ET (δ18O) vs NDVI
row = 2;
[Corr(row,1),Corr(row,2)] = corr(Table7.NDVI, Table7.fT_biased_18O,"Type","Kendall");
lm = fitlm(Table7.NDVI, Table7.fT_biased_18O);
Corr(row,3) = lm.Rsquared.Ordinary;

%T/ET (δ2H) vs LAI
row = 3;
[Corr(row,1),Corr(row,2)] = corr(Table7.LAI, Table7.fT_biased_2H,"Type","Kendall");
lm = fitlm(Table7.LAI, Table7.fT_biased_2H);
Corr(row,3) = lm.Rsquared.Ordinary;

%T/ET (δ2H) vs NDVI
row = 4;
[Corr(row,1),Corr(row,2)] = corr(Table7.NDVI, Table7.fT_biased_2H,"Type","Kendall");
lm = fitlm(Table7.NDVI, Table7.fT_biased_2H);
Corr(row,3) = lm.Rsquared.Ordinary;

Table2D.SiteName = categorical(Table2D.SiteName);
Site_Green.empty = table(1);
clear C ia
[C,ia]           = unique(siteid.site_name);

%Grab each site
for i = 1:length(ia)

Site_Green.(siteid.site_name{ia(i)}) =Table2D(Table2D.SiteName==siteid.site_name{ia(i)},:);

end

%Calculations for each site
for i = 1:length(ia)

%Mean fT by Site
site_mean_fT_18O(i,1) = mean(Site_Green.(siteid.site_name{ia(i)}).fT_biased_18O);
site_mean_fT_2H(i,1)  = mean(Site_Green.(siteid.site_name{ia(i)}).fT_biased_2H);

%Mean NDVI by Site
site_mean_LAI(i,1)  = mean(Site_Green.(siteid.site_name{ia(i)}).LAI);
site_mean_NDVI(i,1) = mean(Site_Green.(siteid.site_name{ia(i)}).NDVI);

end

Table_Site_Green = table(C,site_mean_LAI,site_mean_NDVI,site_mean_fT_18O,site_mean_fT_2H,'VariableNames',{'Site','LAI','NDVI','fT_biased_18O','fT_biased_2H'});

%Corr by Site
%T/ET (δ18O) vs LAI
row = 5;
[Corr(row,1),Corr(row,2)] = corr(Table_Site_Green.LAI, Table_Site_Green.fT_biased_18O,"Type","Kendall");
lm = fitlm(Table_Site_Green.LAI, Table_Site_Green.fT_biased_18O);
Corr(row,3) = lm.Rsquared.Ordinary;

%Inverse
lm = fitlm(Table_Site_Green.LAI, Table_Site_Green.fT_biased_18O,'Weights',1./Table3.fT_biased_SE_18O);
lm = fitlm(Table_Site_Green.NDVI, Table_Site_Green.fT_biased_18O,'Weights',1./Table3.fT_biased_SE_18O);


%T/ET (δ18O) vs NDVI
row = 6;
[Corr(row,1),Corr(row,2)] = corr(Table_Site_Green.NDVI, Table_Site_Green.fT_biased_18O,"Type","Kendall");
lm = fitlm(Table_Site_Green.NDVI, Table_Site_Green.fT_biased_18O);
Corr(row,3) = lm.Rsquared.Ordinary;

%T/ET (δ2H) vs LAI
row = 7;
[Corr(row,1),Corr(row,2)] = corr(Table_Site_Green.LAI, Table_Site_Green.fT_biased_2H,"Type","Kendall");
lm = fitlm(Table_Site_Green.LAI, Table_Site_Green.fT_biased_2H);
Corr(row,3) = lm.Rsquared.Ordinary;

%T/ET (δ2H) vs NDVI
row = 8;
[Corr(row,1),Corr(row,2)] = corr(Table_Site_Green.NDVI, Table_Site_Green.fT_biased_2H,"Type","Kendall");
lm = fitlm(Table_Site_Green.NDVI, Table_Site_Green.fT_biased_2H);
Corr(row,3) = lm.Rsquared.Ordinary;


%% Green & Season 

%Grab each Season

Season_Green.empty = table(1);
Season_Green.Fall  =Table2D(Table2D.Season=="Fall",:);
Season_Green.Summer  =Table2D(Table2D.Season=="Summer",:);

%Inverse Fall
lm = fitlm(Season_Green.Fall.LAI, Season_Green.Fall.fT_biased_18O,'Weights',1./Season_Green.Fall.fT_biased_SE_18O);
lm = fitlm(Season_Green.Fall.NDVI, Season_Green.Fall.fT_biased_18O,'Weights',1./Season_Green.Fall.fT_biased_SE_18O);

%Inverse Summer
lm = fitlm(Season_Green.Summer.LAI, Season_Green.Summer.fT_biased_18O,'Weights',1./Season_Green.Summer.fT_biased_SE_18O);
lm = fitlm(Season_Green.Summer.NDVI, Season_Green.Summer.fT_biased_18O,'Weights',1./Season_Green.Summer.fT_biased_SE_18O);



%%

Table2D.EcosystemType = categorical(Table2D.EcosystemType);

Ecosystem_Green.empty = table(1);

clear D ib
[D,ib] = unique(siteid.ecosystemType);

%Grab each EcosystemType
for i = 1:length(ib)

Ecosystem_Green.(siteid.ecosystemType{ib(i)}) =Table2D(Table2D.EcosystemType==siteid.ecosystemType{ib(i)},:);

end

%Calculations for each ecosystem type
for i = 1:length(ib)

%Mean fT by EcoType
eco_mean_fT_18O(i,1) = mean(Ecosystem_Green.(siteid.ecosystemType{ib(i)}).fT_biased_18O);
eco_mean_fT_2H(i,1)  = mean(Ecosystem_Green.(siteid.ecosystemType{ib(i)}).fT_biased_2H);

%Mean NDVI by Site

eco_mean_LAI(i,1)  = mean(Ecosystem_Green.(siteid.ecosystemType{ib(i)}).LAI);
eco_mean_NDVI(i,1) = mean(Ecosystem_Green.(siteid.ecosystemType{ib(i)}).NDVI);

end

Table_Eco_Green = table(D,eco_mean_LAI,eco_mean_NDVI,eco_mean_fT_18O,eco_mean_fT_2H,'VariableNames',{'Site','LAI','NDVI','fT_biased_18O','fT_biased_2H'});

%Corr by EcosystemType
%T/ET (δ18O) vs LAI
row = 9;
[Corr(row,1),Corr(row,2)] = corr(Table_Eco_Green.LAI, Table_Eco_Green.fT_biased_18O,"Type","Kendall");
lm = fitlm(Table_Eco_Green.LAI, Table_Eco_Green.fT_biased_18O);
Corr(row,3) = lm.Rsquared.Ordinary;

%T/ET (δ18O) vs NDVI
row = 10;
[Corr(row,1),Corr(row,2)] = corr(Table_Eco_Green.NDVI, Table_Eco_Green.fT_biased_18O,"Type","Kendall");
lm = fitlm(Table_Eco_Green.NDVI, Table_Eco_Green.fT_biased_18O);
Corr(row,3) = lm.Rsquared.Ordinary;

%T/ET (δ2H) vs LAI
row = 11;
[Corr(row,1),Corr(row,2)] = corr(Table_Eco_Green.LAI, Table_Eco_Green.fT_biased_2H,"Type","Kendall");
lm = fitlm(Table_Eco_Green.LAI, Table_Eco_Green.fT_biased_2H);
Corr(row,3) = lm.Rsquared.Ordinary;

%T/ET (δ2H) vs NDVI
row = 12;
[Corr(row,1),Corr(row,2)] = corr(Table_Eco_Green.NDVI, Table_Eco_Green.fT_biased_2H,"Type","Kendall");
lm = fitlm(Table_Eco_Green.NDVI, Table_Eco_Green.fT_biased_2H);
Corr(row,3) = lm.Rsqua,red.Ordinary;


%Table 6
Table6 = table(Corr(:,1),Corr(:,2),Corr(:,3),'VariableNames',{'tau','pval','R2'});


%Precipitation from Daymet
P_out    = readtable('P_out.xlsx');

Table7 = addvars(Table7,P_out.PrecipSum,'NewVariableNames','Precip_Ante');

dS_18O     = readtable('dS_18O_summ_out');
dSval_18O  = dS_18O.Var1;
dS_std_18O = (dS_18O.Var2).*sqrt(dS_18O.Var3);

Table8 = table(dET_18O_1day_summ_out.LinSlope,dET_18O_1day_summ_out.LinSE,dET_18O_3day_summ_out.LinSlope,dET_18O_3day_summ_out.LinSE,dET_18O_5day_summ_out.LinSlope,dET_18O_5day_summ_out.LinSE,dET_18O_7day_summ_out.LinSlope,dET_18O_7day_summ_out.LinSE,dETval_comp_18O.mean,dET_SEtotal_18O.std,dSval_18O,dS_std_18O,dEval_18O,dE_SE_18O,'VariableNames',{'dET_1day','dET_1day_SE','dET_3day','dET_3day_SE','dET_5day','dET_5day_SE','dET_7day','dET_7day_SE','dET_Comp','dET_Comp_SE','dS','dS_SE','dE','dE_SE'});

dS_2H     = readtable('dS_2H_summ_out');
dSval_2H  = dS_2H.Var1;
dS_std_2H = (dS_2H.Var2).*sqrt(dS_2H.Var3);

Table9  = table(siteid.site_id,dET_2H_1day_summ_out.LinSlope,dET_2H_1day_summ_out.LinSE,dET_2H_3day_summ_out.LinSlope,dET_2H_3day_summ_out.LinSE,dET_2H_5day_summ_out.LinSlope,dET_2H_5day_summ_out.LinSE,dET_2H_7day_summ_out.LinSlope,dET_2H_7day_summ_out.LinSE,dETval_comp_2H.mean,dET_SEtotal_2H.std,dSval_2H,dS_std_2H,dEval_2H,dE_SE_2H,'VariableNames',{'SiteID','dET_1day','dET_1day_SE','dET_3day','dET_3day_SE','dET_5day','dET_5day_SE','dET_7day','dET_7day_SE','dET_Comp','dET_Comp_SE','dS','dS_SE','dE','dE_SE'});
Table10 = table(siteid.site_id,dETval_comp_18O.mean,dET_SEtotal_18O.std,dEval_18O,dE_SE_18O,dTval_18O,dT_SE_18O,fTval_18O,fT_SE_18O,fTval_18O_biased,fT_SE_18O_biased.std,'VariableNames',{'SiteID','dET_Comp','dET_Comp_SE','dE','dE_SE','dT','dT_SE','fT_NoBias','fT_NoBias_SE','fT_WithBias','fT_WithBias_SE'});
Table11 = table(siteid.site_id,dETval_comp_2H.mean,dET_SEtotal_2H.std,dEval_2H,dE_SE_2H,dTval_2H,dT_SE_2H,fTval_2H,fT_SE_2H,fTval_2H_biased,fT_SE_2H_biased.std,'VariableNames',{'SiteID','dET_Comp','dET_Comp_SE','dE','dE_SE','dT','dT_SE','fT_NoBias','fT_NoBias_SE','fT_WithBias','fT_WithBias_SE'});


Table7a = addvars(Table7,fT_SE_18O_biased.std,fT_SE_2H_biased.std,'NewVariableNames',{'fT_18O_Bias_SE','fT_2H_Bias_SE'});
fT_Precip_Range = nan(9,6);

row = 1;
fT_Precip_Range(row,1) = 0;
fT_Precip_Range(row,2) =numel(Table7a(Table7a.Precip_Ante==0,"fT_biased_18O"));
fT_Precip_Range(row,3) =table2array(mean(Table7a(Table7a.Precip_Ante==0,"fT_biased_18O")));
fT_Precip_Range(row,4) =table2array(sqrt(sum(Table7a(Table7a.Precip_Ante==0,"fT_18O_Bias_SE").^2)./numel(Table7a(Table7a.Precip_Ante==0,"fT_biased_18O"))));
fT_Precip_Range(row,5) =table2array(mean(Table7a(Table7a.Precip_Ante==0,"fT_biased_2H")));
fT_Precip_Range(row,6) =table2array(sqrt(sum(Table7a(Table7a.Precip_Ante==0,"fT_2H_Bias_SE").^2)./numel(Table7a(Table7a.Precip_Ante==0,"fT_biased_2H"))));

row = 2;
fT_Precip_Range(row,1) = 5;
fT_Precip_Range(row,2) =numel(Table7a(Table7a.Precip_Ante<5,"fT_biased_18O"));
fT_Precip_Range(row,3) =table2array(mean(Table7a(Table7a.Precip_Ante<5,"fT_biased_18O")));
fT_Precip_Range(row,4) =table2array(sqrt(sum(Table7a(Table7a.Precip_Ante<5,"fT_18O_Bias_SE").^2)./numel(Table7a(Table7a.Precip_Ante==0,"fT_biased_18O"))));
fT_Precip_Range(row,5) =table2array(mean(Table7a(Table7a.Precip_Ante<5,"fT_biased_2H")));
fT_Precip_Range(row,6) =table2array(sqrt(sum(Table7a(Table7a.Precip_Ante<5,"fT_2H_Bias_SE").^2)./numel(Table7a(Table7a.Precip_Ante==0,"fT_biased_2H"))));

fT_Precip_Range(3:end,1) = [0 5 10 20 30 40 50]';

for i = 3:9

fT_Precip_Range(i,2) =numel(Table7a(Table7a.Precip_Ante>fT_Precip_Range(i,1),"fT_biased_18O"));
fT_Precip_Range(i,3) =table2array(mean(Table7a(Table7a.Precip_Ante>fT_Precip_Range(i,1),"fT_biased_18O")));
fT_Precip_Range(i,4) =table2array(sqrt(sum(Table7a(Table7a.Precip_Ante>fT_Precip_Range(i,1),"fT_18O_Bias_SE").^2)./numel(Table7a(Table7a.Precip_Ante==0,"fT_biased_18O"))));
fT_Precip_Range(i,5) =table2array(mean(Table7a(Table7a.Precip_Ante>fT_Precip_Range(i,1),"fT_biased_2H")));
fT_Precip_Range(i,6) =table2array(sqrt(sum(Table7a(Table7a.Precip_Ante>fT_Precip_Range(i,1),"fT_2H_Bias_SE").^2)./numel(Table7a(Table7a.Precip_Ante==0,"fT_biased_2H"))));

end

Table13 = table(fT_Precip_Range(:,1),fT_Precip_Range(:,2),fT_Precip_Range(:,3),fT_Precip_Range(:,4),fT_Precip_Range(:,5),fT_Precip_Range(:,5),'VariableNames',{'Pre_Thresh','n','fT_18O','fT_18O_SE','fT_2H','fT_2H_SE'});

%Table 14
fT_Precip_Bins = nan(5,6);
fT_Precip_Bins(1,:) = fT_Precip_Range(1,:);

row = 2;
fT_Precip_Bins(row,1) = 0;
fT_Precip_Bins(row,2) =numel(Table7a(Table7a.Precip_Ante<5&Table7a.Precip_Ante>0,"fT_biased_18O"));
fT_Precip_Bins(row,3) =table2array(mean(Table7a(Table7a.Precip_Ante<5&Table7a.Precip_Ante>0,"fT_biased_18O")));
fT_Precip_Bins(row,4) =table2array(sqrt(sum(Table7a(Table7a.Precip_Ante<5&Table7a.Precip_Ante>0,"fT_18O_Bias_SE").^2)./numel(Table7a(Table7a.Precip_Ante==0,"fT_biased_18O"))));
fT_Precip_Bins(row,5) =table2array(mean(Table7a(Table7a.Precip_Ante<5&Table7a.Precip_Ante>0,"fT_biased_2H")));
fT_Precip_Bins(row,6) =table2array(sqrt(sum(Table7a(Table7a.Precip_Ante<5&Table7a.Precip_Ante>0,"fT_2H_Bias_SE").^2)./numel(Table7a(Table7a.Precip_Ante==0,"fT_biased_2H"))));

row = 3;
fT_Precip_Bins(row,1) = 5;
fT_Precip_Bins(row,2) =numel(Table7a(Table7a.Precip_Ante<15&Table7a.Precip_Ante>5,"fT_biased_18O"));
fT_Precip_Bins(row,3) =table2array(mean(Table7a(Table7a.Precip_Ante<15&Table7a.Precip_Ante>5,"fT_biased_18O")));
fT_Precip_Bins(row,4) =table2array(sqrt(sum(Table7a(Table7a.Precip_Ante<15&Table7a.Precip_Ante>5,"fT_18O_Bias_SE").^2)./numel(Table7a(Table7a.Precip_Ante==0,"fT_biased_18O"))));
fT_Precip_Bins(row,5) =table2array(mean(Table7a(Table7a.Precip_Ante<15&Table7a.Precip_Ante>5,"fT_biased_2H")));
fT_Precip_Bins(row,6) =table2array(sqrt(sum(Table7a(Table7a.Precip_Ante<15&Table7a.Precip_Ante>5,"fT_2H_Bias_SE").^2)./numel(Table7a(Table7a.Precip_Ante==0,"fT_biased_2H"))));

row = 4;
fT_Precip_Bins(row,1) = 15;
fT_Precip_Bins(row,2) =numel(Table7a(Table7a.Precip_Ante<35&Table7a.Precip_Ante>15,"fT_biased_18O"));
fT_Precip_Bins(row,3) =table2array(mean(Table7a(Table7a.Precip_Ante<35&Table7a.Precip_Ante>15,"fT_biased_18O")));
fT_Precip_Bins(row,4) =table2array(sqrt(sum(Table7a(Table7a.Precip_Ante<35&Table7a.Precip_Ante>15,"fT_18O_Bias_SE").^2)./numel(Table7a(Table7a.Precip_Ante==0,"fT_biased_18O"))));
fT_Precip_Bins(row,5) =table2array(mean(Table7a(Table7a.Precip_Ante<35&Table7a.Precip_Ante>15,"fT_biased_2H")));
fT_Precip_Bins(row,6) =table2array(sqrt(sum(Table7a(Table7a.Precip_Ante<35&Table7a.Precip_Ante>15,"fT_2H_Bias_SE").^2)./numel(Table7a(Table7a.Precip_Ante==0,"fT_biased_2H"))));


row = 5;
fT_Precip_Bins(row,1) = 35;
fT_Precip_Bins(row,2) =numel(Table7a(Table7a.Precip_Ante>35,"fT_biased_18O"));
fT_Precip_Bins(row,3) =table2array(mean(Table7a(Table7a.Precip_Ante>35,"fT_biased_18O")));
fT_Precip_Bins(row,4) =table2array(sqrt(sum(Table7a(Table7a.Precip_Ante>35,"fT_18O_Bias_SE").^2)./numel(Table7a(Table7a.Precip_Ante>35,"fT_biased_18O"))));
fT_Precip_Bins(row,5) =table2array(mean(Table7a(Table7a.Precip_Ante>35,"fT_biased_2H")));
fT_Precip_Bins(row,6) =table2array(sqrt(sum(Table7a(Table7a.Precip_Ante>35,"fT_2H_Bias_SE").^2)./numel(Table7a(Table7a.Precip_Ante>35,"fT_biased_2H"))));

Table14 = table(fT_Precip_Bins(:,1),fT_Precip_Bins(:,2),fT_Precip_Bins(:,3),fT_Precip_Bins(:,4),fT_Precip_Bins(:,5),fT_Precip_Bins(:,6),'VariableNames',{'Pre_Bins','n','fT_18O','fT_18O_SE','fT_2H','fT_2H_SE'});

Table7b = [Table7a EFlux];


%%

%NIWOT ET MODIS

u = 14;
%ET MODIS Interpolate
clear ET_TT
NEON.NIWO.ET.MOD16A2GF_061_ET_500m
ET_TT = timetable(NEON.(siteid.site_name{u}).ET.Date,(NEON.(siteid.site_name{u}).ET.MOD16A2GF_061_ET_500m),'VariableNames',{'ET_MODIS'});
ET_TT = sortrows(ET_TT);
ET_TT = unique(ET_TT);
NEON.(siteid.site_name{u}).ET_TT_daily = retime(ET_TT,'daily','spline');

clear matching ET_interp_subset
matching = ismember(dateshift(NEON.(siteid.site_name{u}).ET_TT_daily.Time,'start','day'), dateshift(Table2d.SampleDate(ismember(Table2d.SiteName, siteid.site_name{u})),'start','day'));
ET_interp_subset = (NEON.(siteid.site_name{u}).ET_TT_daily.ET_MODIS(matching))/8;

Table7.ET_GapFilled(14:15)  = ET_interp_subset;
Table7a.ET_GapFilled(14:15) = ET_interp_subset;
Table7b.ET_GapFilled(14:15) = ET_interp_subset;

%Adding fT_BIAS to Table10 & 11

Table10 = addvars(Table10,BIAS_fT_18O,'Before','fT_WithBias');
Table11 = addvars(Table11,BIAS_fT_2H,'Before','fT_WithBias');


%% Global Average

mean(Table2.fT_biased_18O)
sqrt(sum(Table2.fT_biased_SE_18O.^2)/56)

mean(Table2.fT_biased_2H)
sqrt(sum(Table2.fT_biased_SE_2H.^2)/56)

% TET and LAI

[tau, pval] = corr(Table7.Precip_Ante,Table7.fT_biased_18O,"Type","Kendall");
lm = fitlm(tbl,'MPG~Weight+Acceleration');

fitlm(Table7,'fT_biased_18O~Precip_Ante+NDVI');

ONAQ_Table7 = Table7(47:56,:);

fitlm(ONAQ_Table7,'fT_biased_18O~Precip_Ante+NDVI')
[tau, pval] = corr(ONAQ_Table7.Precip_Ante,ONAQ_Table7.fT_biased_18O,"Type","Kendall");


% T/ET and ET

[tau, pval] = corr(Table7.ET_Combined,Table7.fT_biased_18O, 'rows','complete',"Type","Kendall");
lm = fitlm(Table7.ET_Combined,Table7.fT_biased_18O);

% Figure 3C&D Inverse Weight

lm = fitlm(Table7b.LAI, Table7b.TFlux_18O,"Weights",1./Table7b.fT_18O_Bias_SE);
lm = fitlm(Table7b.NDVI, Table7b.TFlux_18O,"Weights",1./Table7b.fT_18O_Bias_SE);

%% Saving Tables

%Save output
writetable(Table2,'Tables_Out.xlsx','Sheet','Table2');
writetable(Table2a,'Tables_Out.xlsx','Sheet','Table2a');
writetable(Table2b,'Tables_Out.xlsx','Sheet','Table2b');
writetable(Table2c,'Tables_Out.xlsx','Sheet','Table2c');
writetable(Table2d,'Tables_Out.xlsx','Sheet','Table2d');
writetable(Table2D,'Tables_Out.xlsx','Sheet','Table2D');
writetable(Table3,'Tables_Out.xlsx','Sheet','Table3');
writetable(Table4,'Tables_Out.xlsx','Sheet','Table4');
writetable(Table5,'Tables_Out.xlsx','Sheet','Table5');
writetable(Table6,'Tables_Out.xlsx','Sheet','Table6');
writetable(Table7,'Tables_Out.xlsx','Sheet','Table7');
writetable(Table7a,'Tables_Out.xlsx','Sheet','Table7a');
writetable(Table7b,'Tables_Out.xlsx','Sheet','Table7b');
writetable(Table8,'Tables_Out.xlsx','Sheet','Table8');
writetable(Table9,'Tables_Out.xlsx','Sheet','Table9');

writetable(Table_Eco_Green,'Tables_Out.xlsx','Sheet','Table_Eco_Green');
writetable(Table_Site_Green,'Tables_Out.xlsx','Sheet','Table_Site_Green');


excelFilename = 'SiteLevel.xlsx';
structFieldnames = fieldnames(Site); 
for k = 1:length(structFieldnames)
    fieldname = structFieldnames{k};
    writetable(Site.(fieldname), excelFilename, 'Sheet', sprintf('%s_matlab', fieldname));
end

excelFilename = 'SiteLevelGreen.xlsx';
structFieldnames = fieldnames(Site_Green); 
for k = 1:length(structFieldnames)
    fieldname = structFieldnames{k};
    writetable(Site_Green.(fieldname), excelFilename, 'Sheet', sprintf('%s_matlab', fieldname));
end


excelFilename = 'EcosystemLevel.xlsx';
structFieldnames = fieldnames(Ecosystem); 
for k = 1:length(structFieldnames)
    fieldname = structFieldnames{k};
    writetable(Ecosystem.(fieldname), excelFilename, 'Sheet', sprintf('%s_matlab', fieldname));
end


excelFilename = 'EcosystemLevelGreen.xlsx';
structFieldnames = fieldnames(Ecosystem_Green); 
for k = 1:length(structFieldnames)
    fieldname = structFieldnames{k};
    writetable(Ecosystem_Green.(fieldname), excelFilename, 'Sheet', sprintf('%s_matlab', fieldname));
end

excelFilename = 'Season.xlsx';
structFieldnames = fieldnames(Season); 
for k = 1:length(structFieldnames)
    fieldname = structFieldnames{k};
    writetable(Season.(fieldname), excelFilename, 'Sheet', sprintf('%s_matlab', fieldname));
end

excelFilename = 'SeasonGreen.xlsx';
structFieldnames = fieldnames(Season_Green); 
for k = 1:length(structFieldnames)
    fieldname = structFieldnames{k};
    writetable(Season_Green.(fieldname), excelFilename, 'Sheet', sprintf('%s_matlab', fieldname));
end

excelFilename = 'ONEFLUX.xlsx';
structFieldnames = fieldnames(ONEFLUX); 
for k = 1:length(structFieldnames)
    fieldname = structFieldnames{k};
    writetimetable(ONEFLUX.(fieldname).Flux15_DD_TT, excelFilename, 'Sheet', sprintf('%s_matlab', fieldname));
end


excelFilename = 'LAI_TT.xlsx';
structFieldnames = fieldnames(NEON); 
for k = 1:length(structFieldnames)
    fieldname = structFieldnames{k};
    writetimetable(NEON.(fieldname).LAI_TT_daily, excelFilename, 'Sheet', sprintf('%s_matlab', fieldname));
end

excelFilename = 'NDVI_TT.xlsx';
structFieldnames = fieldnames(NEON); 
for k = 1:length(structFieldnames)
    fieldname = structFieldnames{k};
    writetimetable(NEON.(fieldname).EVI_TT_daily, excelFilename, 'Sheet', sprintf('%s_matlab', fieldname));
end

writetimetable(NEON.NIWO.ET_TT_daily, 'NIWO_ET_MODIS_TT.xlsx', 'Sheet', 'NIWO_matlab');
