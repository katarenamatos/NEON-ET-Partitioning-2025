%% ET Partitioning Code 
%  Code for isotope-based ET partitioning analysis at 13 NEON sites. 
%  Manuscript Title: Evapotranspiration Partitioning Across US Ecoregions: a Multi-Site Study Using Field Stable-Isotope Observations
%  Date last updated: 02/14/2025

%% Gathering Data 
code_path = 'C:\PhD-UNR\Manuscripts\Katarena Matos\PhD Dissertation\Chapter 1\Submission to GRL\Public Code';
cd (code_path)

data_path = 'C:\Users\katar\Documents\Nevada\NEON\NEON-DATA';

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
h5temp = dir ('*.h5')
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
    WaterIso_path = 'C:\Users\katar\Documents\Nevada\NEON\NEON-DATA\NEON_Soil & Plant Water Isotopes_ProjectID 00384\Updated_2021-and-2020';
    cd (fullfile(WaterIso_path))

WIDB_Table              = readtable('1691001063-data.csv');

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
dT_18O_summ_out    = readmatrix('dT_18O_summ_out.xlsx');
dT_2H_summ_out     = readmatrix('dT_2H_summ_out.xlsx');

%% 1-day ET

dET_18O_1day_summ_out = readmatrix('Keeling_Summary_Outputs_tw1day.xlsx','Sheet','18O');
dET_2H_1day_summ_out  = readmatrix('Keeling_Summary_Outputs_tw1day.xlsx','Sheet','2H');

%
D_18O = 0.9727; D_2H  = 0.9757;

for sampdate =  1

%Initiate empty arrays    
Ts = (zeros(1,nn))'; Ta = (zeros(1,nn))'; RH = (zeros(1,nn))'; hA = (zeros(1,nn))'; n = (zeros(1,nn))';

dS_18O = (zeros(1,nn))'; dS_2H = (zeros(1,nn))'; dA_18O = (zeros(1,nn))'; dA_2H = (zeros(1,nn))';
dT_18O = (zeros(1,nn))'; dT_2H = (zeros(1,nn))'; dE_18O = (zeros(1,nn))'; dE_2H = (zeros(1,nn))';
fT_18O = (zeros(1,nn))'; fT_2H = (zeros(1,nn))';


for i = 1:nn

%Sample variable from distribution with mean and std. dev of calculated above   
Ts_temp       = normrnd(Ts_summ_out(sampdate,1),Ts_summ_out(sampdate,2)); Ts(i) = Ts_temp;
Ta_temp       = normrnd(Ta_summ_out(sampdate,1),Ta_summ_out(sampdate,2)); Ta(i) = Ta_temp;
RH_temp       = normrnd(RH_summ_out(sampdate,1),RH_summ_out(sampdate,2)); RH(i) = RH_temp;

dS_18O_temp   = normrnd(dSoil_18O_summ_out(sampdate,1),dSoil_18O_summ_out(sampdate,2)); dS_18O (i) = dS_18O_temp;
dS_2H_temp    = normrnd(dSoil_2H_summ_out(sampdate,1),dSoil_2H_summ_out(sampdate,2));   dS_2H (i) = dS_2H_temp;
dA_18O_temp   = normrnd(dA_18O_summ_out(sampdate,1),dA_18O_summ_out(sampdate,2));       dA_18O (i) = dA_18O_temp;
dA_2H_temp    = normrnd(dA_2H_summ_out(sampdate,1),dA_2H_summ_out(sampdate,2));         dA_2H (i) = dA_2H_temp;

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

% Variable (6): n 
%Sample from normal distribution over possible range of "n"
n_temp = 0.5 + (1-0.5)*rand(1,1); n(i) = n_temp;

% Variable (7): εk
ek_18O_temp = n_temp*(1 - hA_temp)*(1-D_18O);
ek_2H_temp  = n_temp*(1 - hA_temp)*(1-D_2H);

%Craig and Gordon (1965)
dE_18O_temp = (((alpha_18O_temp*dS_18O_temp-hA_temp*dA_18O_temp))-(eps_eq_18O_temp+ek_18O_temp)/((1-hA_temp)+(ek_18O_temp/1000))); dE_18O (i) = dE_18O_temp;
dE_2H_temp  = (((alpha_2H_temp*dS_2H_temp-hA_temp*dA_2H_temp))-(eps_eq_2H_temp+ek_18O_temp)/((1-hA_temp)+(ek_2H_temp/1000)));  dE_2H (i) = dE_2H_temp;

%Sample variable from distribution with mean and std. dev of calculated above   
dT_18O_temp = normrnd(dT_18O_summ_out(sampdate,1),dT_18O_summ_out(sampdate,2)); dT_18O (i) = dT_18O_temp;
dT_2H_temp  = normrnd(dT_2H_summ_out(sampdate,1),dT_2H_summ_out(sampdate,2)); dT_2H (i) = dT_2H_temp;

dET_18O_temp = dET_18O_1day_summ_out(sampdate,2); dET_18O (i) = dET_18O_temp;
dET_2H_temp  = dET_2H_1day_summ_out(sampdate,2); dET_2H (i) = dET_2H_temp;

end

dET_18O = dET_18O';
dET_2H  = dET_2H';

% Calculate T/ET 
fT_18O_temp = (dET_18O-dE_18O./(dT_18O-dE_18O)); fT_18O = fT_18O_temp;
fT_2H_temp  = (dET_2H-dE_2H)./(dT_2H-dE_2H);     fT_2H  = fT_2H_temp;

output_concat = horzcat(Ts, Ta, RH, hA, n, dS_18O, dS_2H, dA_18O, dA_2H, dT_18O, dT_2H, dE_18O, dE_2H, fT_18O, fT_2H);
output_concat_table = array2table(output_concat,"VariableNames",{'Ts', 'Ta', 'RH', 'hA', 'n', 'dS_18O', 'dS_2H', 'dA_18O', 'dA_2H', 'dT_18O', 'dT_2H', 'dE_18O', 'dE_2H', 'fT_18O', 'fT_2H'});
writetable(output_concat_table,'All-Vars_MC-Output_Matlab_dET_3day.xlsx','Sheet',siteid.site_id{sampdate});

fTET_18O_summout(sampdate,1) = mean(fT_18O);
fTET_18O_summout(sampdate,2) = std(fT_18O);
fTET_18O_summout(sampdate,3) = nn;

fTET_2H_summout(sampdate,1) = mean(fT_2H);
fTET_2H_summout(sampdate,2) = std(fT_2H);
fTET_2H_summout(sampdate,3) = nn;

end

cd (code_path)
writematrix(fTET_18O_summout,'MC_Output_Matlab_nn1000_1dET_fTET_18O_summout.xlsx','Sheet',[siteid.site_id{sampdate} char(' fT_18O')]);
writematrix(fTET_2H_summout,'MC_Output_Matlab_nn1000_1dET_fTET_2H_summout.xlsx','Sheet',[siteid.site_id{sampdate} char(' fT_2H')]);


%% 3-day ET

dET_18O_3day_summ_out = readmatrix('Keeling_Summary_Outputs_tw3day.xlsx','Sheet','18O');
dET_2H_3day_summ_out  = readmatrix('Keeling_Summary_Outputs_tw3day.xlsx','Sheet','2H');

%
D_18O = 0.9727; D_2H  = 0.9757;

for sampdate =  1

%Initiate empty arrays    
Ts = (zeros(1,nn))'; Ta = (zeros(1,nn))'; RH = (zeros(1,nn))'; hA = (zeros(1,nn))'; n = (zeros(1,nn))';

dS_18O = (zeros(1,nn))'; dS_2H = (zeros(1,nn))'; dA_18O = (zeros(1,nn))'; dA_2H = (zeros(1,nn))';
dT_18O = (zeros(1,nn))'; dT_2H = (zeros(1,nn))'; dE_18O = (zeros(1,nn))'; dE_2H = (zeros(1,nn))';
fT_18O = (zeros(1,nn))'; fT_2H = (zeros(1,nn))';


for i = 1:nn

%Sample variable from distribution with mean and std. dev of calculated above   
Ts_temp       = normrnd(Ts_summ_out(sampdate,1),Ts_summ_out(sampdate,2)); Ts(i) = Ts_temp;
Ta_temp       = normrnd(Ta_summ_out(sampdate,1),Ta_summ_out(sampdate,2)); Ta(i) = Ta_temp;
RH_temp       = normrnd(RH_summ_out(sampdate,1),RH_summ_out(sampdate,2)); RH(i) = RH_temp;

dS_18O_temp   = normrnd(dSoil_18O_summ_out(sampdate,1),dSoil_18O_summ_out(sampdate,2)); dS_18O (i) = dS_18O_temp;
dS_2H_temp    = normrnd(dSoil_2H_summ_out(sampdate,1),dSoil_2H_summ_out(sampdate,2));   dS_2H (i) = dS_2H_temp;
dA_18O_temp   = normrnd(dA_18O_summ_out(sampdate,1),dA_18O_summ_out(sampdate,2));       dA_18O (i) = dA_18O_temp;
dA_2H_temp    = normrnd(dA_2H_summ_out(sampdate,1),dA_2H_summ_out(sampdate,2));         dA_2H (i) = dA_2H_temp;

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

% Variable (6): n 
%Sample from normal distribution over possible range of "n"
n_temp = 0.5 + (1-0.5)*rand(1,1); n(i) = n_temp;

% Variable (7): εk
ek_18O_temp = n_temp*(1 - hA_temp)*(1-D_18O);
ek_2H_temp  = n_temp*(1 - hA_temp)*(1-D_2H);

%Craig and Gordon (1965)
dE_18O_temp = (((alpha_18O_temp*dS_18O_temp-hA_temp*dA_18O_temp))-(eps_eq_18O_temp+ek_18O_temp)/((1-hA_temp)+(ek_18O_temp/1000))); dE_18O (i) = dE_18O_temp;
dE_2H_temp  = (((alpha_2H_temp*dS_2H_temp-hA_temp*dA_2H_temp))-(eps_eq_2H_temp+ek_18O_temp)/((1-hA_temp)+(ek_2H_temp/1000)));  dE_2H (i) = dE_2H_temp;

%Sample variable from distribution with mean and std. dev of calculated above   
dT_18O_temp = normrnd(dT_18O_summ_out(sampdate,1),dT_18O_summ_out(sampdate,2)); dT_18O (i) = dT_18O_temp;
dT_2H_temp  = normrnd(dT_2H_summ_out(sampdate,1),dT_2H_summ_out(sampdate,2)); dT_2H (i) = dT_2H_temp;

dET_18O_temp = dET_18O_3day_summ_out(sampdate,2); dET_18O (i) = dET_18O_temp;
dET_2H_temp  = dET_2H_3day_summ_out(sampdate,2); dET_2H (i) = dET_2H_temp;

end

dET_18O = dET_18O';
dET_2H  = dET_2H';

% Calculate T/ET 
fT_18O_temp = (dET_18O-dE_18O./(dT_18O-dE_18O)); fT_18O = fT_18O_temp;
fT_2H_temp  = (dET_2H-dE_2H)./(dT_2H-dE_2H);     fT_2H  = fT_2H_temp;

output_concat = horzcat(Ts, Ta, RH, hA, n, dS_18O, dS_2H, dA_18O, dA_2H, dT_18O, dT_2H, dE_18O, dE_2H, fT_18O, fT_2H);
output_concat_table = array2table(output_concat,"VariableNames",{'Ts', 'Ta', 'RH', 'hA', 'n', 'dS_18O', 'dS_2H', 'dA_18O', 'dA_2H', 'dT_18O', 'dT_2H', 'dE_18O', 'dE_2H', 'fT_18O', 'fT_2H'});
writetable(output_concat_table,'All-Vars_MC-Output_Matlab_dET_3day.xlsx','Sheet',siteid.site_id{sampdate});


fTET_18O_summout(sampdate,1) = mean(fT_18O);
fTET_18O_summout(sampdate,2) = std(fT_18O);
fTET_18O_summout(sampdate,3) = nn;

fTET_2H_summout(sampdate,1) = mean(fT_2H);
fTET_2H_summout(sampdate,2) = std(fT_2H);
fTET_2H_summout(sampdate,3) = nn;

end

cd (code_path)
writematrix(fTET_18O_summout,'MC_Output_Matlab_nn1000_3dET_fTET_18O_summout.xlsx','Sheet',[siteid.site_id{sampdate} char(' fT_18O')]);
writematrix(fTET_2H_summout,'MC_Output_Matlab_nn1000_3dET_fTET_2H_summout.xlsx','Sheet',[siteid.site_id{sampdate} char(' fT_2H')]);


%% 5-day ET

dET_18O_5day_summ_out = readmatrix('Keeling_Summary_Outputs_tw5day.xlsx','Sheet','18O');
dET_2H_5day_summ_out  = readmatrix('Keeling_Summary_Outputs_tw5day.xlsx','Sheet','2H');

%
D_18O = 0.9727; D_2H  = 0.9757;

for sampdate =  1

%Initiate empty arrays    
Ts = (zeros(1,nn))'; Ta = (zeros(1,nn))'; RH = (zeros(1,nn))'; hA = (zeros(1,nn))'; n = (zeros(1,nn))';

dS_18O = (zeros(1,nn))'; dS_2H = (zeros(1,nn))'; dA_18O = (zeros(1,nn))'; dA_2H = (zeros(1,nn))';
dT_18O = (zeros(1,nn))'; dT_2H = (zeros(1,nn))'; dE_18O = (zeros(1,nn))'; dE_2H = (zeros(1,nn))';
fT_18O = (zeros(1,nn))'; fT_2H = (zeros(1,nn))';


for i = 1:nn

%Sample variable from distribution with mean and std. dev of calculated above   
Ts_temp       = normrnd(Ts_summ_out(sampdate,1),Ts_summ_out(sampdate,2)); Ts(i) = Ts_temp;
Ta_temp       = normrnd(Ta_summ_out(sampdate,1),Ta_summ_out(sampdate,2)); Ta(i) = Ta_temp;
RH_temp       = normrnd(RH_summ_out(sampdate,1),RH_summ_out(sampdate,2)); RH(i) = RH_temp;

dS_18O_temp   = normrnd(dSoil_18O_summ_out(sampdate,1),dSoil_18O_summ_out(sampdate,2)); dS_18O (i) = dS_18O_temp;
dS_2H_temp    = normrnd(dSoil_2H_summ_out(sampdate,1),dSoil_2H_summ_out(sampdate,2));   dS_2H (i) = dS_2H_temp;
dA_18O_temp   = normrnd(dA_18O_summ_out(sampdate,1),dA_18O_summ_out(sampdate,2));       dA_18O (i) = dA_18O_temp;
dA_2H_temp    = normrnd(dA_2H_summ_out(sampdate,1),dA_2H_summ_out(sampdate,2));         dA_2H (i) = dA_2H_temp;

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

% Variable (6): n 
%Sample from normal distribution over possible range of "n"
n_temp = 0.5 + (1-0.5)*rand(1,1); n(i) = n_temp;

% Variable (7): εk
ek_18O_temp = n_temp*(1 - hA_temp)*(1-D_18O);
ek_2H_temp  = n_temp*(1 - hA_temp)*(1-D_2H);

%Craig and Gordon (1965)
dE_18O_temp = (((alpha_18O_temp*dS_18O_temp-hA_temp*dA_18O_temp))-(eps_eq_18O_temp+ek_18O_temp)/((1-hA_temp)+(ek_18O_temp/1000))); dE_18O (i) = dE_18O_temp;
dE_2H_temp  = (((alpha_2H_temp*dS_2H_temp-hA_temp*dA_2H_temp))-(eps_eq_2H_temp+ek_18O_temp)/((1-hA_temp)+(ek_2H_temp/1000)));  dE_2H (i) = dE_2H_temp;

%Sample variable from distribution with mean and std. dev of calculated above   
dT_18O_temp = normrnd(dT_18O_summ_out(sampdate,1),dT_18O_summ_out(sampdate,2)); dT_18O (i) = dT_18O_temp;
dT_2H_temp  = normrnd(dT_2H_summ_out(sampdate,1),dT_2H_summ_out(sampdate,2)); dT_2H (i) = dT_2H_temp;

dET_18O_temp = dET_18O_5day_summ_out(sampdate,2); dET_18O (i) = dET_18O_temp;
dET_2H_temp  = dET_2H_5day_summ_out(sampdate,2); dET_2H (i) = dET_2H_temp;

end

dET_18O = dET_18O';
dET_2H  = dET_2H';

% Calculate T/ET 
fT_18O_temp = (dET_18O-dE_18O./(dT_18O-dE_18O)); fT_18O = fT_18O_temp;
fT_2H_temp  = (dET_2H-dE_2H)./(dT_2H-dE_2H);     fT_2H  = fT_2H_temp;

output_concat = horzcat(Ts, Ta, RH, hA, n, dS_18O, dS_2H, dA_18O, dA_2H, dT_18O, dT_2H, dE_18O, dE_2H, fT_18O, fT_2H);
output_concat_table = array2table(output_concat,"VariableNames",{'Ts', 'Ta', 'RH', 'hA', 'n', 'dS_18O', 'dS_2H', 'dA_18O', 'dA_2H', 'dT_18O', 'dT_2H', 'dE_18O', 'dE_2H', 'fT_18O', 'fT_2H'});
writetable(output_concat_table,'All-Vars_MC-Output_Matlab_dET_5day.xlsx','Sheet',siteid.site_id{sampdate});

fTET_18O_summout(sampdate,1) = mean(fT_18O);
fTET_18O_summout(sampdate,2) = std(fT_18O);
fTET_18O_summout(sampdate,3) = nn;

fTET_2H_summout(sampdate,1) = mean(fT_2H);
fTET_2H_summout(sampdate,2) = std(fT_2H);
fTET_2H_summout(sampdate,3) = nn;

end

cd (code_path)
writematrix(fTET_18O_summout,'MC_Output_Matlab_nn1000_5dET_fTET_18O_summout.xlsx','Sheet',[siteid.site_id{sampdate} char(' fT_18O')]);
writematrix(fTET_2H_summout,'MC_Output_Matlab_nn1000_5dET_fTET_2H_summout.xlsx','Sheet',[siteid.site_id{sampdate} char(' fT_2H')]);




%% 7-day ET

dET_18O_7day_summ_out = readmatrix('Keeling_Summary_Outputs_tw7day.xlsx','Sheet','18O');
dET_2H_7day_summ_out  = readmatrix('Keeling_Summary_Outputs_tw7day.xlsx','Sheet','2H');

%
D_18O = 0.9727; D_2H  = 0.9757;

for sampdate =  1

%Initiate empty arrays    
Ts = (zeros(1,nn))'; Ta = (zeros(1,nn))'; RH = (zeros(1,nn))'; hA = (zeros(1,nn))'; n = (zeros(1,nn))';

dS_18O = (zeros(1,nn))'; dS_2H = (zeros(1,nn))'; dA_18O = (zeros(1,nn))'; dA_2H = (zeros(1,nn))';
dT_18O = (zeros(1,nn))'; dT_2H = (zeros(1,nn))'; dE_18O = (zeros(1,nn))'; dE_2H = (zeros(1,nn))';
fT_18O = (zeros(1,nn))'; fT_2H = (zeros(1,nn))';


for i = 1:nn

%Sample variable from distribution with mean and std. dev of calculated above   
Ts_temp       = normrnd(Ts_summ_out(sampdate,1),Ts_summ_out(sampdate,2)); Ts(i) = Ts_temp;
Ta_temp       = normrnd(Ta_summ_out(sampdate,1),Ta_summ_out(sampdate,2)); Ta(i) = Ta_temp;
RH_temp       = normrnd(RH_summ_out(sampdate,1),RH_summ_out(sampdate,2)); RH(i) = RH_temp;

dS_18O_temp   = normrnd(dSoil_18O_summ_out(sampdate,1),dSoil_18O_summ_out(sampdate,2)); dS_18O (i) = dS_18O_temp;
dS_2H_temp    = normrnd(dSoil_2H_summ_out(sampdate,1),dSoil_2H_summ_out(sampdate,2));   dS_2H (i) = dS_2H_temp;
dA_18O_temp   = normrnd(dA_18O_summ_out(sampdate,1),dA_18O_summ_out(sampdate,2));       dA_18O (i) = dA_18O_temp;
dA_2H_temp    = normrnd(dA_2H_summ_out(sampdate,1),dA_2H_summ_out(sampdate,2));         dA_2H (i) = dA_2H_temp;

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

% Variable (6): n 
%Sample from normal distribution over possible range of "n"
n_temp = 0.5 + (1-0.5)*rand(1,1); n(i) = n_temp;

% Variable (7): εk
ek_18O_temp = n_temp*(1 - hA_temp)*(1-D_18O);
ek_2H_temp  = n_temp*(1 - hA_temp)*(1-D_2H);

%Craig and Gordon (1965)
dE_18O_temp = (((alpha_18O_temp*dS_18O_temp-hA_temp*dA_18O_temp))-(eps_eq_18O_temp+ek_18O_temp)/((1-hA_temp)+(ek_18O_temp/1000))); dE_18O (i) = dE_18O_temp;
dE_2H_temp  = (((alpha_2H_temp*dS_2H_temp-hA_temp*dA_2H_temp))-(eps_eq_2H_temp+ek_18O_temp)/((1-hA_temp)+(ek_2H_temp/1000)));  dE_2H (i) = dE_2H_temp;

%Sample variable from distribution with mean and std. dev of calculated above   
dT_18O_temp = normrnd(dT_18O_summ_out(sampdate,1),dT_18O_summ_out(sampdate,2)); dT_18O (i) = dT_18O_temp;
dT_2H_temp  = normrnd(dT_2H_summ_out(sampdate,1),dT_2H_summ_out(sampdate,2)); dT_2H (i) = dT_2H_temp;

dET_18O_temp = dET_18O_7day_summ_out(sampdate,2); dET_18O (i) = dET_18O_temp;
dET_2H_temp  = dET_2H_7day_summ_out(sampdate,2); dET_2H (i) = dET_2H_temp;

end

dET_18O = dET_18O';
dET_2H  = dET_2H';

% Calculate T/ET 
fT_18O_temp = (dET_18O-dE_18O./(dT_18O-dE_18O)); fT_18O = fT_18O_temp;
fT_2H_temp  = (dET_2H-dE_2H)./(dT_2H-dE_2H);     fT_2H  = fT_2H_temp;

output_concat = horzcat(Ts, Ta, RH, hA, n, dS_18O, dS_2H, dA_18O, dA_2H, dT_18O, dT_2H, dE_18O, dE_2H, fT_18O, fT_2H);
output_concat_table = array2table(output_concat,"VariableNames",{'Ts', 'Ta', 'RH', 'hA', 'n', 'dS_18O', 'dS_2H', 'dA_18O', 'dA_2H', 'dT_18O', 'dT_2H', 'dE_18O', 'dE_2H', 'fT_18O', 'fT_2H'});
writetable(output_concat_table,'All-Vars_MC-Output_Matlab_dET_7day.xlsx','Sheet',siteid.site_id{sampdate});

fTET_18O_summout(sampdate,1) = mean(fT_18O);
fTET_18O_summout(sampdate,2) = std(fT_18O);
fTET_18O_summout(sampdate,3) = nn;

fTET_2H_summout(sampdate,1) = mean(fT_2H);
fTET_2H_summout(sampdate,2) = std(fT_2H);
fTET_2H_summout(sampdate,3) = nn;

end

cd (code_path)
writematrix(fTET_18O_summout,'MC_Output_Matlab_nn1000_7dET_fTET_18O_summout.xlsx','Sheet',[siteid.site_id{sampdate} char(' fT_18O')]);
writematrix(fTET_2H_summout,'MC_Output_Matlab_nn1000_7dET_fTET_2H_summout.xlsx','Sheet',[siteid.site_id{sampdate} char(' fT_2H')]);


