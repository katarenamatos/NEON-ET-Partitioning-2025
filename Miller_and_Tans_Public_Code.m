%% ET Partitioning Code: Keeling-Type Analysis for dET in 1-,3-,5-, and 7- time windows
%  Code for isotope-based ET partitioning analysis at 13 NEON sites. 
%  Manuscript Title: Evapotranspiration Partitioning Across US Ecoregions: a Multi-Site Study Using Field Stable-Isotope Observations
%  Date last updated: 02/14/2025

% Gathering Data 
code_path = 'C:\PhD-UNR\Manuscripts\Katarena Matos\PhD Dissertation\Chapter 1\Submission to GRL\Public Code';
cd (code_path)

%Load NEON sites Metadata Including ID, Sample Date, and NEON Domains 
siteid      = readtable("siteid.xlsx");
SampDate    = readtable("SampDate_start_and_end_day.xlsx");
site_date   = readtable("site_date.xlsx");
site_domain = readmatrix("site_domain.xlsx");

%% dET 
% Description: Isotopic composition of dET from isotopic composition of atmosphere 

% Data Product used below: NEON H2O Calibrated 
    % Set path to location of this data product on the local computer

data_path = 'C:\Users\katar\Documents\NEON-Local\NEON H2O Calibrated';
cd (data_path)

%NEON file naming convention
filname1    = 'NEON.D%02.f.%s*';

%Eddy Covariance Levels
level = {'Lvl1','Lvl2','Lvl3','Lvl4','Lvl5','Lvl6','Lvl7','Lvl8'};

%Extract correct data file for each site and sample date
u = 1 

site = string(siteid.site_name{u});
filname2    = sprintf (filname1, site_domain(u), siteid.site_name{u});
filname     = dir(fullfile(pwd,filname2));
info_h5 = h5info(filname.name)

clear s 
s.Lvl0 = timetable(datetime('today'),1); %initiate structure

for l = 2:eddy_levels(u) %Loop through levels of Eddy Towers 

path_iso_9m     = getfield(info_h5.Groups.Groups.Groups(1).Groups.Groups,{l},'Name');

%Iso 18O
clear data_iso18O_9m data_iso2H_9m data_rtioMoleWet_9m date stor_LvL
path_iso_9m_18O = append(path_iso_9m,'/dlta18OH2o');
data_iso18O_9m  = h5read(filname.name,path_iso_9m_18O);

%Iso 2H
path_iso_9m_2H      = append(path_iso_9m,'/dlta2HH2o');
data_iso2H_9m = h5read(filname.name,path_iso_9m_2H);

%Ratio Mole Wet H2O
path_iso_9m_rtioMoleWetH2o      = append(path_iso_9m,'/rtioMoleWetH2o');
data_rtioMoleWet_9m = h5read(filname.name,path_iso_9m_rtioMoleWetH2o);

%Manipulate dates
data_iso18O_9m.timeBgn = data_iso18O_9m.timeBgn';
data_iso18O_9m.timeEnd = data_iso18O_9m.timeEnd';

for i =1:length(data_iso18O_9m.timeEnd)
    data_iso18O_9m.timeEndNEW(i,1) = datetime(join(data_iso18O_9m.timeEnd(i,1:end),''),'InputFormat','yyyy-MM-dd''T''HH:mm:ss.sss''Z');
    data_iso18O_9m.timeBgnNEW(i,1) = datetime(join(data_iso18O_9m.timeBgn(i,1:end),''),'InputFormat','yyyy-MM-dd''T''HH:mm:ss.sss''Z');
end

date = data_iso18O_9m.timeEndNEW;
stor_LvL = timetable(date,data_iso18O_9m.mean,data_iso2H_9m.mean,data_rtioMoleWet_9m.mean,'VariableNames',{'dlta18OH2o','dlta2HH2o','rtioMoleWetH2o'});

s.(level{l}) = stor_LvL;
str = append(site,' ',level(l))

%cd (code_path)
%writetimetable(stor_LvL,'NEON_Compiled_Keeling_Inputs.xlsx','Sheet',str);

end

%Save outputs to structure
S.(site).hdf5 = timetable(datetime('today'),1); 
S.(site) = s;

%% Keeling Type: 1-day Time Window 

%time window
tw = {'tw_1day'} %,'3day','5day','7day'}

%dummer varible
w = 1;

%Loop through sample dates, compute dET
for u = 1 %:size(SampDate,1)

site = string(siteid.site_name{u});
clear X Y_18O Y_2H tr 

%1-day
tr    = timerange(SampDate.Start_Sample_Day(u),SampDate.End_Sample_Day(u));
X     = nan(30,eddy_levels(u));
Y_18O = nan(30,eddy_levels(u));
Y_2H  = nan(30,eddy_levels(u));

clear X Y2H Y18O test2 test3
for l = 1:eddy_levels(u) %Loop through levels of Eddy Towers 

clear test x y18O y2H 
%for each Eddy Covriance Level, extract sample date of interest
test = S.(site).(level{l})(tr,:);

%Keeling-Type Regression Variables
x    = test.rtioMoleWetH2o;
y18O = test.dlta18OH2o.*x;
y2H  = test.dlta2HH2o.*x;

X(1:length(x),l)    = x;
Y18O(1:length(x),l) = y18O;
Y2H(1:length(x),l)  = y2H;

%save regression variables into an array
%for a given sample date, all levels are stacked/concatenated into 'test2'
if l == 1
test2 = table2array(test);
else
test2 = vertcat(test2, table2array(test));
end

end 

%Columns 1, 2, and 3 correspond to: dlta18OH2O; dlta2HH2O, rtioMoleWetH2O 

test3 = test2;

%Get rid of any rows with NaN 
clear x y_18O y_2H fitlm_18O fitlm_2H p18O p2H yy_18O yy2H
test3(any(isnan(test3), 2), :) = [];
x     = test3(:,3);
y_18O = test3(:,1).*x;
y_2H  = test3(:,2).*x;

%Linear regression of Keeping-Type variables (18O)
fitlm_18O = fitlm(x,y_18O);
p18O      = polyfit(x,y_18O,1);
x_eval    = linspace(0,max(x)+5,40);
yy_18O    = polyval(p18O,x_eval);

%Save stats for each regression (18O)
stats_18O_tw1day(u,1) = fitlm_18O.Rsquared.Ordinary;
stats_18O_tw1day (u,2) = fitlm_18O.Coefficients.pValue(2);
stats_18O_tw1day (u,3) = fitlm_18O.Coefficients.Estimate(2); %linear slope
stats_18O_tw1day (u,4) = fitlm_18O.Coefficients.SE(2);

%Linear regression of Keeping-Type variables (2H)
fitlm_2H = fitlm(x,y_2H);
p2H      = polyfit(x,y_2H,1);
yy_2H    = polyval(p2H,x_eval);
x_eval   = linspace(0,max(x)+5,40);
yy_2H    = polyval(p2H,x_eval);

%Save stats for each regression (2H)
stats_2H_tw1day (u,1) = fitlm_2H.Rsquared.Ordinary;
stats_2H_tw1day (u,2) = fitlm_2H.Coefficients.pValue(2);
stats_2H_tw1day (u,3) = fitlm_2H.Coefficients.Estimate(2); %linear slope
stats_2H_tw1day (u,4) = fitlm_2H.Coefficients.SE(2);

%Save Keeling-type variables/data used in regression
output_table = array2table(test3,"VariableNames",{'dlta18OH2O','dlta2HH2O','rtioMoleWetH2O'});

%Save outputs
cd (code_path)
writetable(output_table,'Compiled_Tower.xlsx','Sheet',[siteid.site_id{u},'+',char(tw{w})]);

end 

%Save Outputs
cd (code_path)
output_stats18O_tw1day =array2table(stats_18O_tw1day,'VariableNames',{'Lin R2','pval','Lin Slope','Lin SE'});
output_stats18O_tw1day =addvars(output_stats18O_tw1day,siteid.site_id,'Before','Lin R2','NewVariableNames','site_id');
writetable(output_stats18O_tw1day,'Keeling_Summary_Outputs.xlsx','Sheet','18O')

output_stats2H_tw1day =array2table(stats_2H_tw1day,'VariableNames',{'Lin R2','pval','Lin Slope','Lin SE'});
output_stats2H_tw1day =addvars(output_stats2H_tw1day,siteid.site_id,'Before','Lin R2','NewVariableNames','site_id');
writetable(output_stats2H_tw1day,'Keeling_Summary_Outputs.xlsx','Sheet','2H')
writetable(output_stats2H_tw1day,'Keeling_Summary_Outputs.xlsx','Sheet','2H')
%% Keeling Type: 3-day Time Window 

%time window
tw = {'3day'};

%dummer varible
w = 1;

%Loop through sample dates, compute dET
for u = 1 %:size(SampDate,1)

site = string(siteid.site_name{u});
clear X Y_18O Y_2H tr 

%1-day
tr    = timerange(SampDate.Start_Sample_Day(u)-1,SampDate.End_Sample_Day(u)+1);
X     = nan(30,eddy_levels(u));
Y_18O = nan(30,eddy_levels(u));
Y_2H  = nan(30,eddy_levels(u));

clear X Y2H Y18O test2 test3
for l = 1:eddy_levels(u) %Loop through levels of Eddy Towers 

clear test x y18O y2H 
%for each Eddy Covriance Level, extract sample date of interest
test = S.(site).(level{l})(tr,:);

%Keeling-Type Regression Variables
x    = test.rtioMoleWetH2o;
y18O = test.dlta18OH2o.*x;
y2H  = test.dlta2HH2o.*x;

X(1:length(x),l)    = x;
Y18O(1:length(x),l) = y18O;
Y2H(1:length(x),l)  = y2H;

%save regression variables into an array
%for a given sample date, all levels are stacked/concatenated into 'test2'
if l == 1
test2 = table2array(test);
else
test2 = vertcat(test2, table2array(test));
end

end 

%Columns 1, 2, and 3 correspond to: dlta18OH2O; dlta2HH2O, rtioMoleWetH2O 

test3 = test2;

%Get rid of any rows with NaN 
clear x y_18O y_2H fitlm_18O fitlm_2H p18O p2H yy_18O yy2H
test3(any(isnan(test3), 2), :) = [];
x     = test3(:,3);
y_18O = test3(:,1).*x;
y_2H  = test3(:,2).*x;

%Linear regression of Keeping-Type variables (18O)
fitlm_18O = fitlm(x,y_18O);
p18O      = polyfit(x,y_18O,1);
x_eval    = linspace(0,max(x)+5,40);
yy_18O    = polyval(p18O,x_eval);

%Save stats for each regression (18O)
stats_18O_tw3day (u,1) = fitlm_18O.Rsquared.Ordinary;
stats_18O_tw3day (u,2) = fitlm_18O.Coefficients.pValue(2);
stats_18O_tw3day (u,3) = fitlm_18O.Coefficients.Estimate(2); %linear slope
stats_18O_tw3day (u,4) = fitlm_18O.Coefficients.SE(2);

%Linear regression of Keeping-Type variables (2H)
fitlm_2H = fitlm(x,y_2H);
p2H      = polyfit(x,y_2H,1);
yy_2H    = polyval(p2H,x_eval);
x_eval   = linspace(0,max(x)+5,40);
yy_2H    = polyval(p2H,x_eval);

%Save stats for each regression (2H)
stats_2H_tw3day (u,1) = fitlm_2H.Rsquared.Ordinary;
stats_2H_tw3day (u,2) = fitlm_2H.Coefficients.pValue(2);
stats_2H_tw3day (u,3) = fitlm_2H.Coefficients.Estimate(2); %linear slope
stats_2H_tw3day (u,4) = fitlm_2H.Coefficients.SE(2);

%Save Keeling-type variables/data used in regression
output_table = array2table(test3,"VariableNames",{'dlta18OH2O','dlta2HH2O','rtioMoleWetH2O'});

%Save outputs
cd (code_path)
writetable(output_table,'Compiled_Tower.xlsx','Sheet',[siteid.site_id{u},'+',char(tw{w})]);

end 

%Save Outputs
cd (code_path)
output_stats18O_tw3day =array2table(stats_18O_tw3day,'VariableNames',{'Lin R2','pval','Lin Slope','Lin SE'});
output_stats18O_tw3day =addvars(output_stats18O_tw3day,siteid.site_id,'Before','Lin R2','NewVariableNames','site_id');
writetable(output_stats18O_tw3day,'Keeling_Summary_Outputs.xlsx','Sheet','18O')

output_stats2H_tw3day =array2table(stats_2H_tw3day,'VariableNames',{'Lin R2','pval','Lin Slope','Lin SE'});
output_stats2H_tw3day =addvars(output_stats2H_tw3day,siteid.site_id,'Before','Lin R2','NewVariableNames','site_id');
writetable(output_stats2H_tw3day,'Keeling_Summary_Outputs.xlsx','Sheet','2H')
writetable(output_stats2H_tw3day,'Keeling_Summary_Outputs.xlsx','Sheet','2H')

%% Keeling Type: 5-day Time Window 

%time window
tw = {'5day'};

%dummer varible
w = 1;

%Loop through sample dates, compute dET
for u = 1 %:size(SampDate,1)

site = string(siteid.site_name{u});
clear X Y_18O Y_2H tr 

%1-day
tr    = timerange(SampDate.Start_Sample_Day(u)-2,SampDate.End_Sample_Day(u)+2);
X     = nan(30,eddy_levels(u));
Y_18O = nan(30,eddy_levels(u));
Y_2H  = nan(30,eddy_levels(u));

clear X Y2H Y18O test2 test3
for l = 1:eddy_levels(u) %Loop through levels of Eddy Towers 

clear test x y18O y2H 
%for each Eddy Covriance Level, extract sample date of interest
test = S.(site).(level{l})(tr,:);

%Keeling-Type Regression Variables
x    = test.rtioMoleWetH2o;
y18O = test.dlta18OH2o.*x;
y2H  = test.dlta2HH2o.*x;

X(1:length(x),l)    = x;
Y18O(1:length(x),l) = y18O;
Y2H(1:length(x),l)  = y2H;

%save regression variables into an array
%for a given sample date, all levels are stacked/concatenated into 'test2'
if l == 1
test2 = table2array(test);
else
test2 = vertcat(test2, table2array(test));
end

end 

%Columns 1, 2, and 3 correspond to: dlta18OH2O; dlta2HH2O, rtioMoleWetH2O 

test3 = test2;

%Get rid of any rows with NaN 
clear x y_18O y_2H fitlm_18O fitlm_2H p18O p2H yy_18O yy2H
test3(any(isnan(test3), 2), :) = [];
x     = test3(:,3);
y_18O = test3(:,1).*x;
y_2H  = test3(:,2).*x;

%Linear regression of Keeping-Type variables (18O)
fitlm_18O = fitlm(x,y_18O);
p18O      = polyfit(x,y_18O,1);
x_eval    = linspace(0,max(x)+5,40);
yy_18O    = polyval(p18O,x_eval);

%Save stats for each regression (18O)
stats_18O_tw5day (u,1) = fitlm_18O.Rsquared.Ordinary;
stats_18O_tw5day (u,2) = fitlm_18O.Coefficients.pValue(2);
stats_18O_tw5day (u,3) = fitlm_18O.Coefficients.Estimate(2); %linear slope
stats_18O_tw5day (u,4) = fitlm_18O.Coefficients.SE(2);

%Linear regression of Keeping-Type variables (2H)
fitlm_2H = fitlm(x,y_2H);
p2H      = polyfit(x,y_2H,1);
yy_2H    = polyval(p2H,x_eval);
x_eval   = linspace(0,max(x)+5,40);
yy_2H    = polyval(p2H,x_eval);

%Save stats for each regression (2H)
stats_2H_tw5day (u,1) = fitlm_2H.Rsquared.Ordinary;
stats_2H_tw5day (u,2) = fitlm_2H.Coefficients.pValue(2);
stats_2H_tw5day (u,3) = fitlm_2H.Coefficients.Estimate(2); %linear slope
stats_2H_tw5day (u,4) = fitlm_2H.Coefficients.SE(2);

%Save Keeling-type variables/data used in regression
output_table = array2table(test3,"VariableNames",{'dlta18OH2O','dlta2HH2O','rtioMoleWetH2O'});

%Save outputs
cd (code_path)
writetable(output_table,'Compiled_Tower.xlsx','Sheet',[siteid.site_id{u},'+',char(tw{w})]);

end 

cd (code_path)
output_stats18O_tw5day =array2table(stats_18O_tw5day,'VariableNames',{'Lin R2','pval','Lin Slope','Lin SE'});
output_stats18O_tw5day =addvars(output_stats18O_tw5day,siteid.site_id,'Before','Lin R2','NewVariableNames','site_id');
writetable(output_stats18O_tw5day,'Keeling_Summary_Outputs.xlsx','Sheet','18O')

output_stats2H_tw5day =array2table(stats_2H_tw5day,'VariableNames',{'Lin R2','pval','Lin Slope','Lin SE'});
output_stats2H_tw5day =addvars(output_stats2H_tw5day,siteid.site_id,'Before','Lin R2','NewVariableNames','site_id');
writetable(output_stats2H_tw5day,'Keeling_Summary_Outputs.xlsx','Sheet','2H')
writetable(output_stats2H_tw5day,'Keeling_Summary_Outputs.xlsx','Sheet','2H')

%% Keeling Type: 7-day Time Window 

%time window
tw = {'7day'};

%dummer varible
w = 1;

%Loop through sample dates, compute dET
for u = 1 %:size(SampDate,1)

site = string(siteid.site_name{u});
clear X Y_18O Y_2H tr 

%1-day
tr    = timerange(SampDate.Start_Sample_Day(u)-3,SampDate.End_Sample_Day(u)+3);
X     = nan(30,eddy_levels(u));
Y_18O = nan(30,eddy_levels(u));
Y_2H  = nan(30,eddy_levels(u));

clear X Y2H Y18O test2 test3
for l = 1:eddy_levels(u) %Loop through levels of Eddy Towers 

clear test x y18O y2H 
%for each Eddy Covriance Level, extract sample date of interest
test = S.(site).(level{l})(tr,:);

%Keeling-Type Regression Variables
x    = test.rtioMoleWetH2o;
y18O = test.dlta18OH2o.*x;
y2H  = test.dlta2HH2o.*x;

X(1:length(x),l)    = x;
Y18O(1:length(x),l) = y18O;
Y2H(1:length(x),l)  = y2H;

%save regression variables into an array
%for a given sample date, all levels are stacked/concatenated into 'test2'
if l == 1
test2 = table2array(test);
else
test2 = vertcat(test2, table2array(test));
end

end 

%Columns 1, 2, and 3 correspond to: dlta18OH2O; dlta2HH2O, rtioMoleWetH2O 

test3 = test2;

%Get rid of any rows with NaN 
clear x y_18O y_2H fitlm_18O fitlm_2H p18O p2H yy_18O yy2H
test3(any(isnan(test3), 2), :) = [];
x     = test3(:,3);
y_18O = test3(:,1).*x;
y_2H  = test3(:,2).*x;

%Linear regression of Keeping-Type variables (18O)
fitlm_18O = fitlm(x,y_18O);
p18O      = polyfit(x,y_18O,1);
x_eval    = linspace(0,max(x)+5,40);
yy_18O    = polyval(p18O,x_eval);

%Save stats for each regression (18O)
stats_18O_tw7day (u,1) = fitlm_18O.Rsquared.Ordinary;
stats_18O_tw7day (u,2) = fitlm_18O.Coefficients.pValue(2);
stats_18O_tw7day (u,3) = fitlm_18O.Coefficients.Estimate(2); %linear slope
stats_18O_tw7day (u,4) = fitlm_18O.Coefficients.SE(2);

%Linear regression of Keeping-Type variables (2H)
fitlm_2H = fitlm(x,y_2H);
p2H      = polyfit(x,y_2H,1);
yy_2H    = polyval(p2H,x_eval);
x_eval   = linspace(0,max(x)+5,40);
yy_2H    = polyval(p2H,x_eval);

%Save stats for each regression (2H)
stats_2H_tw7day (u,1) = fitlm_2H.Rsquared.Ordinary;
stats_2H_tw7day (u,2) = fitlm_2H.Coefficients.pValue(2);
stats_2H_tw7day (u,3) = fitlm_2H.Coefficients.Estimate(2); %linear slope
stats_2H_tw7day (u,4) = fitlm_2H.Coefficients.SE(2);

%Save Keeling-type variables/data used in regression
output_table = array2table(test3,"VariableNames",{'dlta18OH2O','dlta2HH2O','rtioMoleWetH2O'});

%Save outputs
cd (code_path)
writetable(output_table,'Compiled_Tower.xlsx','Sheet',[siteid.site_id{u},'+',char(tw{w})]);

end 

%Save Outputs
cd (code_path)
output_stats18O_tw7day =array2table(stats_18O_tw7day,'VariableNames',{'Lin R2','pval','Lin Slope','Lin SE'});
output_stats18O_tw7day =addvars(output_stats18O_tw7day,siteid.site_id(1),'Before','Lin R2','NewVariableNames','site_id');
writetable(output_stats18O_tw7day,'Keeling_Summary_Outputs_tw7day.xlsx','Sheet','18O')

output_stats2H_tw7day =array2table(stats_2H_tw7day,'VariableNames',{'Lin R2','pval','Lin Slope','Lin SE'});
output_stats2H_tw7day =addvars(output_stats2H_tw7day,siteid.site_id,'Before','Lin R2','NewVariableNames','site_id');
writetable(output_stats2H_tw7day,'Keeling_Summary_Outputs_tw7day.xlsx','Sheet','2H')
