%%
% fitler_sensitivty_exp_analysis
subjID = {'004';'006';'009';'011';'100'};%'011',

pathResults = '/Users/jossando/trabajo/image_filtering/06_RawData';

acuities = [];
gaussian_dprimes = [];
gaussian_01cutoffhalf_dprimes = [];
gaussian_01cutoff_dprimes = [];
butterworth_custom_01cutoff_order3_dprimes = [];
butterworth_custom_01cutoff_order6_dprimes = [];
CSF_filter_dprimes = [];

for ss = 1:length(subjID)
    load(fullfile(pathResults,['FS' subjID{ss}]))
    acuities = [acuities;PARAMS.subject.visualAcuity];
    gaussian_dprimes = [gaussian_dprimes;TRIALS.gaussian_dprime];
    gaussian_01cutoffhalf_dprimes = [gaussian_01cutoffhalf_dprimes;TRIALS.gaussian_01cutoffhalf_dprime];
    gaussian_01cutoff_dprimes = [gaussian_01cutoff_dprimes;TRIALS.gaussian_01cutoff_dprime];
    butterworth_custom_01cutoff_order3_dprimes = [butterworth_custom_01cutoff_order3_dprimes;TRIALS.butterworth_custom_01cutoff_order3_dprime];
    butterworth_custom_01cutoff_order6_dprimes = [butterworth_custom_01cutoff_order6_dprimes;TRIALS.butterworth_custom_01cutoff_order6_dprime];
    CSF_filter_dprimes = [CSF_filter_dprimes;TRIALS.CSF_filter_dprime];
end
resultTable = table(subjID,acuities,gaussian_dprimes, gaussian_01cutoffhalf_dprimes,gaussian_01cutoff_dprimes,butterworth_custom_01cutoff_order3_dprimes,butterworth_custom_01cutoff_order6_dprimes,CSF_filter_dprimes);
4
filtNames = {'gauss','gauss01half','guass01','butt0103','butt0106','csf'};
figjp,hold on
for ss = 1:height(resultTable)
    scatter(1:6,table2array(resultTable(ss,3:end)),50,'o')
    plot(1:6,table2array(resultTable(ss,3:end)),'Color',[.9 .9 .9],'LineWidth',.5)
end
axis([0 7 -1 4])
set(gca,'xTick',1:6,'XTickLabels',filtNames)
ylabel('d prime')
hline(0)
