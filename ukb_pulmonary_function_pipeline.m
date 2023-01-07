%% SNPs Associated w/ Lung Function Among AFR and EUR

%% Defination of Variables 
% *Forced vital capacity (FVC):* is the amount of air that can be forcibly
% exhaled from your lungs after taking the deepest breath possible, as
% measured by spirometry. This test may help distinguish obstructive lung
% diseases, such as asthma and COPD, from restrictive lung diseases, such
% as pulmonary fibrosis and sarcoidosis.

% *Forced Expiratory Volume (FEV1):* is the maximum amount of air you can
% forcefully exhale in one second. It is used to describe the degree of
% airway obstruction caused by asthma in a routine test called spirometry
% or pulmonary function testing, using an instrument called a spirometer.

% *Peak Expiratory Flow Rate (PEFR):* is the maximum flow rate generated
% during a forceful exhalation, starting from full lung inflation.  PEFR
% primarily reflects large airway flow and depends on the voluntary effort
% and muscular strength of the patient.

clc; clear; close all force

try
    cd '/Users/sinkala/Documents/MATLAB/UKBioBank/Manuscript - Respiratory'
catch
    fprintf('\nWorking on cluster or a different computer\n')
end

% add the path to the GWAS results for Africans and Europeans 
addpath(['/Users/sinkala/Documents/MATLAB/UKBioBank/' ...
    'Manuscript - Respiratory/eur_afr_gwas'])

% download the gwas summary statistics from AWS at the following links 
% FCV: https://pan-ukb-us-east-1.s3.amazonaws.com/sumstats_flat_files/continuous-3062-both_sexes-irnt.tsv.bgz 
% FEV1: https://pan-ukb-us-east-1.s3.amazonaws.com/sumstats_flat_files/continuous-3063-both_sexes-irnt.tsv.bgz 
% PEF: https://pan-ukb-us-east-1.s3.amazonaws.com/sumstats_flat_files/continuous-3064-both_sexes-irnt.tsv.bgz

% unzip then untar the AWS files and rename them as follows: 
% ForcedExpiratoryVolume_gwas.txt
% ForcedVitalCapacity_gwas.txt and 
% PeakExpiratoryFlow_gwas.txt

%% Load the UK Biobank data

try
    fprintf('\n Now loading the data \n')
    data = readtable("hpc_ukbiobank_clinical.csv");
    data.EthnicGroupBW = categorical(data.EthnicGroupBW);
catch
    error( sprintf(['\nYou will need to download the clinical ' ...
        'information files from the UK biobank:\n\nThese file include:\n' ...
        '1. ukb40800.csv\n' ...
        '2. Data_Dictionary_Showcase.csv\n' ...
        '3. Codings_Showcase\n' ...
        '\nUpon downloading the following file run the script called \n'...
        'ukBioBankDataProcessing.m which will produce the ', ...
        'hpc_ukbiobank_clinical.csv file required to continue \nthe ' ...
        'analysis using this code\n']) )
end

% get smaller chuck of the data
data = data(:,{'EthnicGroupBW','ForcedVitalCapacity_FVC_', ...
    'ForcedExpiratoryVolumeIn1_second_FEV1_','PeakExpiratoryFlow_PEF_',...
    'BodyMassIndex_BMI_','AgeAtRecruitment','StandingHeight',...
    'WheezeOrWhistlingInTheChestInLastYear',...
    'ShortnessOfBreathWalkingOnLevelGround',...
    'ChestPainOrDiscomfortWalkingNormally',...
    'ChestPainFeltOutsidePhysicalActivity',...
    'ChestPainOrDiscomfortWhenWalkingUphillOrHurrying',...
    'BloodClot_DVT_Bronchitis_Emphysema_Asthma_Rhinitis_Eczema_Aller'} );

%% ***************** Smoking Question by Reviewer ******************

% *********************** First set of comments **********************
% Reviewer 2: I didn't see any mention of smoking in this manuscript
% despite the well-known effects of smoking in pulmonary function measures.
% Failure to account for an effect of smoking could reduce statistical
% power.

% Reviewer 3: Similarly, lung function without adjustment for smoking
% status in GWAS might result in signals associated to smoking behaviour.

% ******************** Second set of comments *************************

% Reviewer 2: The existence of a lot of smoking variables is certainly not
% a reason that not to adjust for smoking. Neither is the fact that it
% might be challenging. By asking for adjustment for smoking, no one was
% suggesting that analyses be adjusted for all 34 smoking variables.
% Analyses such as these are generally adjusted for either current smoking
% status (Yes/No) or separately for both current smoking and ever smoking
% status (yes/no). These variables would be easy to generate from the first
% two smoking variables listed. If the authors are not going to adjust for
% smoking, then a compelling scientific or analytical justification for not
% doing so is necessary. The reasons listed are insufficient.

% Reviwer 3: I agree there could be a lot of work done to investigate the
% impact of smoking behaviour. However, the author could at least adjust
% for current smoking status like many other lung function GWAS did.

% Specify the model using a formula that allows up to two-way interactions
% between the variables age, weight, and sex. Smoker is the response
% variable. ~EthnicGroupBW','+AgeAtRecruitment

% % load the data datasets 
% glmData = readtable('cleaned_UK_Biobank_Blacks_and_Whites.csv');
% 
% myVarsLong = {'ForcedExpiratoryVolumeIn1_second_FEV1_',...
%     'ForcedVitalCapacity_FVC_',...
%     'PeakExpiratoryFlow_PEF_'} ;
% 
% % convert the first to Ethinicity 
% glmData.Properties.VariableNames(1) = "Ethnicity";
% 
% % get only a subset of the data
% glmData = glmData(:,  ['Ethnicity', myVarsLong, ...
%     glmData.Properties.VariableNames( contains( ...
%     glmData.Properties.VariableNames,'smok','IgnoreCase',true) ) ] ) ;
% 
% % run the generalised model hares
% for ii = 1:length(myVarsLong)
% 
%     % here are the model specs 
%     modelAge = sprintf('%s ~ Ethnicity*AgeAtRecruitment',myVarsLong{ii}) ;
% 
%     % Fit a generalised linear model for Age at recruitment
%     glmModel.(shortTitles{ii}).age =  ...
%         fitglm(glmData,modelAge,'Distribution','normal') ; 
% end
% 
% clear myVarsLong modelAge modelHeight ii glmData

%% Process the gwas results from the uk biobank into two set

% one whould be for africans and the other from europeans
ukGwasFiles = {'ForcedExpiratoryVolume_gwas.txt';...
    'ForcedVitalCapacity_gwas.txt';'PeakExpiratoryFlow_gwas.txt'};

% check if the variant qc materic file exist
if ~exist('variant_qc_metrics.txt','file')
    variantqc = readtable('full_variant_qc_metrics.txt') ;
    variantqc = variantqc(:,{'chrom','pos','ref','alt',...
        'rsid','nearest_genes'}) ;
    writetable(variantqc,'variant_qc_metrics.txt') ;
end

% check if the file exists
if ~exist(replace(ukGwasFiles{1},'.txt','_AFR.txt'),'file')
    % load the qc matrics of the data and match the first variable name to
    % that of the gwas results
    variantqc = readtable('variant_qc_metrics.txt') ;
    variantqc.Properties.VariableNames([1,5]) = {'chr','ID'};
    variantqc.alt = [];

    % process the files
    for ii = 1:length(ukGwasFiles)

        % load the gewas results
        curResults = readtable(ukGwasFiles{ii}) ;

        % print something to the screen
        fprintf('\n Creating GWAS summary results for %s: number %d\n',...
            extractBefore(ukGwasFiles{ii},'_gwas'), ii)

        % break it into two files
        gwasAFR = removevars(curResults,'pval_EUR');
        gwasEUR = removevars(curResults,'pval_AFR');

        % change the variable names
        gwasAFR.Properties.VariableNames(end) = "P" ;
        gwasEUR.Properties.VariableNames(end) = "P" ;

        % add the variant types to the table: Variant fields#

        % chr: Chromosome of the variant. pos: Position of the variant in
        % GRCh37 coordinates. ref: Reference allele on the forward strand.
        % alt: Alternate allele (not necessarily minor allele). Used as
        % effect allele for GWAS.
        gwasAFR = innerjoin(gwasAFR,variantqc,'Keys',{'chr','pos','ref'});
        gwasEUR = innerjoin(gwasEUR,variantqc,'Keys',{'chr','pos','ref'});

        % change the variable names to upper case to match those required
        % for plotting
        gwasAFR.Properties.VariableNames = upper( ...
            gwasAFR.Properties.VariableNames) ;
        gwasEUR.Properties.VariableNames = upper( ...
            gwasEUR.Properties.VariableNames) ;

        % throw in an assertion 
        assert( height(gwasAFR) == height(gwasEUR))

        % save to two text files
        writetable(gwasAFR,replace(ukGwasFiles{ii},'.txt','_AFR.txt') );
        writetable(gwasEUR,replace(ukGwasFiles{ii},'.txt','_EUR.txt') );

    end

    % remove the variantqc data which is very large
    clear variantqc

end

%% Run the T-tests: comparing between AFR and EUR

% the colors to use the for two groups
groupColors = [0.00,0.45,0.74; 0.64,0.08,0.18] ;
 
% get the variables
myVars = {'ForcedVitalCapacity_FVC_', ...
    'ForcedExpiratoryVolumeIn1_second_FEV1_',...
    'PeakExpiratoryFlow_PEF_'};

% the name of the plot and the table variable names
figTitles = {'Forced Vital Capacity', ...
    'Forced Expiratory Volume in 1 Second','Peak Expiratory Flow'};

% set up the unit and the short titles of the variable names
varUnits = {'L','L','L/min'} ;
shortTitles = {'FVC','FEV1','PEF'} ;

% preallocated the t-test table
ttestTable = table('Size',[1,7],'VariableTypes',...
    {'string','double','double','double','double','double','double'}, ...
    'VariableNames',{'Variable','meanAFR','meanEUR','tValue', ...
    'lowerBound','upperBound','pValue'});

% run the test in a loop
for ii = 1:length( myVars)
    
    fprintf('\n Running t-test and plotting data for %s \n', ...
        figTitles{ii}) 
    
    % get the current var name
    curVarName = myVars{ii} ;
    
    % get the data in the current column and fill the outliers using linear
    % interpolation
    curData = filloutliers(data.(curVarName),'linear') ;

    % perform a ttest with unequal variable assumed
    [~,p,ci,stats] = ttest2(curData(data.EthnicGroupBW =='AFR'),...
        curData(data.EthnicGroupBW == 'EUR'),'Vartype','unequal') ;

    % get the mean values for the two groups
    meanAFR = mean(curData(data.EthnicGroupBW == 'AFR'),'omitnan');
    meanEUR = mean(curData(data.EthnicGroupBW == 'EUR'),'omitnan');

    % add to the table
    ttestTable.Variable(ii) = figTitles(ii) ;
    ttestTable.meanAFR(ii) = meanAFR;
    ttestTable.meanEUR(ii) = meanEUR ;
    ttestTable.lowerBound(ii) = ci(1);
    ttestTable.upperBound(ii) = ci(2);
    ttestTable.pValue(ii) = p ;
    ttestTable.tValue(ii) = stats.tstat;
    
    % plot the figure
    colourBoxPlot(curData,data.EthnicGroupBW,flipud(groupColors),true);
    hold on 
    
    % put the p value on the plot depending on its value
    if p < 0.0001
        text(0.5, 0.9 , ['p ', convertPValue2SuperScript(p)],...
            'Units','normalized','FontWeight','bold','FontSize',12, ...
            'HorizontalAlignment','center')
    else
        text(0.5, 0.9 ,sprintf('p = %0.4f',p),...
            'Units','normalized','FontWeight','bold','FontSize',12, ...
            'HorizontalAlignment','center')
    end
    
    % add the title to the plot and a y-axis label
    if ii == 3 % for peak expiratory follow
        ylabel([ shortTitles{ii} ,' (L/min)'])
    else
        ylabel([ shortTitles{ii} ,' (L)'])
    end
    title(shortTitles{ii},'FontSize',16,'FontWeight','bold')
    hold off
    
    % save the figures 
    saveas(gcf,[shortTitles{ii},'_ttest.fig'],'fig')

end 

hold off 

clear aa ans ci curData curVarName groups ii meanAFR meanEUR ...
    p plotNames stats xValues plotName 

%% Save the Source Data for Figure 1

% here are the data 
dataS = data(:, 1:7) ;
dataS.Properties.VariableNames(1) = "EthnicGroup" ;
writetable(dataS,'FigureSoureData.xlsx','Sheet','Fig1 Source Data');

clear dataS

%% ******************** REVIEWER COMMENT *******************************

% Please discuss how the 5-year difference in age might have affected these
% findings. 5 years could be significant in terms of lung health – younger
% AFR individuals might mean less variability in the lung function
% measurements, which could have affected power to detect associations as
% well.
[~,p,ci,stats] = ttest2( ...
    data.AgeAtRecruitment(data.EthnicGroupBW =='AFR'),...
    data.AgeAtRecruitment(data.EthnicGroupBW == 'EUR'),...
    'Vartype','unequal') ;

% plot the figure
colourBoxPlot(data.AgeAtRecruitment, data.EthnicGroupBW, ...
    flipud(groupColors),true);
hold on

% put the p value on the plot depending on its value
if p < 0.0001
    text(0.5, 0.9 , ['p ', convertPValue2SuperScript(p)],...
        'Units','normalized','FontWeight','bold','FontSize',12, ...
        'HorizontalAlignment','center')
else
    text(0.5, 0.9 ,sprintf('p = %0.4f',p),...
        'Units','normalized','FontWeight','bold','FontSize',12, ...
        'HorizontalAlignment','center')
end

% add the title to the plot and a y-axis label
ylabel('Age (years)')
title('Age','FontSize',16,'FontWeight','bold')
hold off

% save the figures
saveas(gcf,'Age t-test.fig','fig')

%% Make the Manhattan Plot for the GWAS data 

% come up with the new variable for the gwas analysis 
% get the variables
myVars = {'ForcedVitalCapacity','ForcedExpiratoryVolume',...
    'PeakExpiratoryFlow'};

% here are the snps to exclude

% Name:	SNPs recommended for exclusion from analyses Size:	1,805 MD5:
% cee103826b9ffbacd9d119180bf98249 This document lists a number of
% genotyped autosomal SNPs (65) which have been found to show significantly
% different allele frequencies between the UK BiLEVE array and the UK
% Biobank array. These SNPs are in the interim data release but should be
% excluded from analyses.
% 
% A number (27) of these SNPs were used in phasing and imputation. UKB
% strongly recommends conditioning on array in association tests to
% ameliorate the effect of these SNPs. There could still be a subtle bias
% in the neighbourhood of these SNPs after conditioning, but this will
% depend upon the phenotype being tested for association. UKB recommends
% looking carefully at any results with imputed SNPs in the regions of the
% affected SNPs, including confirming any GWAS hits with the genotyped-only
% data and looking at cluster plots of the genotype data.
% 
% Additionally, there are a number of SNPs (46) on chromosome X which show
% a significant allele frequency difference between males and females or
% show differences between arrays. UKB recommends that these SNPs be
% excluded from all analyses.

    
% for the FEV in One second
for ii = 1:length(myVars)

    fprintf('\n Plotting the Manhattan plot for %s \n', figTitles{ii})

    figure()
    % plot the manhattam plots in a loop
    ManhattanPlot( [myVars{ii}, '_gwas_AFR.txt'], ...
        'title',[shortTitles{ii},' - AFR'])

    figure()
    ManhattanPlot([myVars{ii}, '_gwas_EUR.txt'],...
        'title',[shortTitles{ii} ,' - EUR'])

end

clear sigSnpsAFR sigSnpsEUR 

%% Reviewer Comment

% It is critical that a effect size plot (I would have an estimated beta in
% recent European ancestry vs estimated beta in recent African ancestry -
% both with error bars - plot - this scatter plot will get at whether this
% is mainly an effect size shift - which I doubt or a allele frequency and
% sample size issue - which probably is the case). Remember that the
% estimate of betas is feasible as long as you have at least some
% individuals (it is just that the error bar might be very wide)

% Get the SNPs and the beta values from the PanUKBiobank repository and
% plot them 


%% Find the Common Variant in  EUR and AFR for Each Feature

% Also make some Q-Q plots

% preallocate a table for the significant snps
snpSigTable = table('Size',[1,5],'VariableTypes',...
    {'string','double','double','double','cellstr'}, ...
    'VariableNames',{'Variable','sigSnpsAFR','sigSnpsEUR', ...
    'commonSigSnps','listSigSnps'});

% for the each variable
if ~exist([shortTitles{1},'_venn_all_snps.fig'],'file')

    % run this in a for loop
    for ii = 1:length(myVars)

        fprintf('\nFinding the common SNPs plot for %s \n', figTitles{ii})

        % get the data in a loop
        curDataAFR = readtable([myVars{ii}, '_gwas_AFR.txt']) ;
        curDataEUR = readtable([myVars{ii}, '_gwas_EUR.txt']) ;

        % make the Q-Q plot Create figure
        figure()
        hold on
        pd = makedist('Lognormal');
        subplot(1,2,1)
        qqplot(-log10(curDataAFR.P), pd) % chi square values
        % set some figure properties and add title ot the figure
        set(gca,'FontSize',12,'LineWidth',1,'Box','off')
        ylabel(['Quantiles of ', shortTitles{ii} ])
        title(['Q-Q Plot of ', shortTitles{ii} ,' AFR'],'FontSize',16,...
            'FontWeight','bold')

        subplot(1,2,2)
        qqplot(-log10(curDataEUR.P),pd)
        % set some figure properties and add title ot the figure
        set(gca,'FontSize',12,'LineWidth',1,'Box','off')
        ylabel(['Quantiles of ',shortTitles{ii}])
        title(['Q-Q Plot of ', shortTitles{ii},' EUR' ],'FontSize',16,...
            'FontWeight','bold')
        hold off

        % save the figure
        saveas(gcf,[shortTitles{ii},'_qqplot.fig'],'fig')
        saveas(gcf,[shortTitles{ii},'_qqplot.png'],'png')

        % get only the significant snps from both tables
        curDataAFR = curDataAFR(curDataAFR.P <= 5e-8,:) ;
        curDataEUR = curDataEUR(curDataEUR.P <= 5e-8,:) ;

        % how many snps are significant in EUR and in black for the current
        % variable
        snpSigTable.Variable(ii) = figTitles(ii);
        snpSigTable.sigSnpsAFR(ii) = height(curDataAFR);
        snpSigTable.sigSnpsEUR(ii) = height(curDataEUR) ;

        % get the common snps
        commonSnps = intersect(curDataAFR.ID, curDataEUR.ID) ;
        snpSigTable.commonSigSnps(ii) = length(commonSnps);
        snpSigTable.listSigSnps{ii} = commonSnps' ;

        % draw the venn diagram for the SNPs that are intersecting
        figure()
        venn([200,200],75,'EdgeColor','black') ;
        hold on

        % add text to the figure: get numbers of intersects to add to the
        % plots
        myNumText(1) = height(curDataAFR) - length(commonSnps) ;
        myNumText(2) = length(commonSnps);
        myNumText(3) = height(curDataEUR) - length(commonSnps) ;

        % here are text positions and the titles of the plots
        textPos = [0.25, 0.47 , 0.70] ;
        for jj = 1:length(textPos)
            text(textPos(jj),0.5,num2str(myNumText(jj)),...
                'Units','normalized','FontWeight','bold','FontSize',16,...
                'HorizontalAlignment','center')
        end

        % add the titles to the circle
        myGroups = {'AFR','EUR'};
        textPos = [0.33 , 0.66] ;
        for jj = 1:length(textPos)
            text(textPos(jj), 0.98 ,[shortTitles{ii}, ' ', myGroups{jj}],...
                'Units','normalized','FontWeight','bold',...
                'FontSize',14,'HorizontalAlignment','center')
        end
        hold off

        % save the venn diagram
        saveas(gcf,[shortTitles{ii},'_venn_all_snps.fig'],'fig')

        % sort the rows of the table
        curDataAFR = sortrows(curDataAFR,'P','ascend');
        curDataEUR = sortrows(curDataEUR,'P','ascend');

        % move the variable in the table to make them usable for
        % finemapping
        curDataAFR = movevars(movevars(curDataAFR,'ID','Before',1),'P',...
            'After','POS') ;
        curDataEUR = movevars(movevars(curDataEUR,'ID','Before',1),'P',...
            'After','POS') ;

        % change the p-value
        curDataAFR.P = -log10(curDataAFR.P);
        curDataEUR.P = -log10(curDataEUR.P);

        % save to text files to be use for getting the casual variants
        writetable(curDataAFR(:,1:4),...
            ['sigSnp',shortTitles{ii},'_AFR.txt'], ...
            'delimiter', ' ','WriteVariableNames',false)
        writetable(curDataEUR(:,1:4),...
            ['sigSnp',shortTitles{ii},'_EUR.txt'], ...
            'delimiter', ' ','WriteVariableNames', false')

    end

    disp(snpSigTable)

end

clear ii curDataEUR curDataAFR commonSnps ii  textPos myNumText

%% ************************** Review 2 Comment **************************

% Perhaps a plot comparing the p-values in each ancestry group for all
% variants that were above a certain threshold of statistical significance?

% check if the files with snps that are significant in any group are
% present
if ~exist([myVars{1},'_lung_anySig.txt'],'file')

    % run this in a for loop
    for ii = 1:length(myVars)

        fprintf('\nFinding the common SNPs plot for %s \n', myVars{ii})

        % time the running of this programme
        tstart = tic ;

        % make a data store the data
        dsAFR = tabularTextDatastore([myVars{ii},'_gwas_AFR.txt']);
        dsEUR = tabularTextDatastore([myVars{ii},'_gwas_EUR.txt']);

        % load the snps
        curDataAFR = tall(dsAFR) ;
        curDataEUR = tall(dsEUR) ;

        % get the significant snps in african
        sigAFR = gather(curDataAFR.P <= 5e-8);
        sigEUR = gather(curDataEUR.P <= 5e-8);

        % get only the significant snps from both tables
        curDataAFR = curDataAFR(sigAFR | sigEUR ,:);
        curDataEUR = curDataEUR(sigAFR | sigEUR ,:);

        % change the variable names for the pvalue
        curDataAFR.Properties.VariableNames(5) = "AFR" ;
        curDataEUR.Properties.VariableNames(5) = "EUR" ;

        % get both snps in one table
        bothData = innerjoin(curDataAFR, curDataEUR(:,{'EUR','ID'}), ...
            'Keys','ID');

        % gather the datasets
        bothData = gather(bothData);

        % write the data to a file
        writetable(bothData,[myVars{ii},'_lung_anySig.txt'])

        telapsed = toc(tstart)/60 ;
        fprintf('\n Time elapsed is %s Minutes \n', num2str(telapsed) )
    end
end

% find the snps above a threshold on 0.05 in each groups when they are
% significnt in one group

% preallocate a table and the struct of results 
sigSums = table('Size',[3,5],...
    'VariableTypes',{'string','double','double','double','double'} , ...
    'VariableNames',{'Measure','SigAFR','P<0.05inEUR','SigEUR',...
    'P<0.05inAFR'} ) ;

% run this in a for loop
for ii = 1:length(myVars)

    % read in the data 
    bothData = readtable([myVars{ii},'_lung_anySig.txt']) ;
    bothData = movevars(bothData,'AFR','After','EUR') ;

    % how many snps are above the 0.05 in afr or eur when sig in one group
    sigSums.Measure(ii) = shortTitles(ii) ;
    sigSums.SigAFR(ii) = nnz(bothData.AFR <= 5e-8) ;
    sigSums.SigEUR(ii) = nnz(bothData.EUR <= 5e-8) ;
    sigSums.('P<0.05inEUR')(ii) = nnz(bothData.AFR <= 5e-8 & ...
        bothData.EUR <= 0.05) ;
    sigSums.('P<0.05inAFR')(ii) = nnz(bothData.EUR <= 5e-8 & ...
        bothData.AFR <= 0.05) ;

    % add the data to the structured array 
    bothData = sortrows(bothData,"AFR","ascend");
    sigSumStruct.(shortTitles{ii}) = bothData ;
 
%     % here are the groups 
%     theGroups = {'EUR','AFR'};
% 
%     % here is the plot of the p-values
%     for jj = 1:length(theGroups)
% 
%         % get the current ploting data
%         curPlotData = bothData(:,{'CHR','POS',theGroups{jj}}) ;
%         curPlotData.Properties.VariableNames(2:end) = ["BP","P"] ;
% 
%         % plot the manhattan figures
%         ManhattanPlotGrey(curPlotData, 'plotColor', groupColors(jj,:), ...
%             'title','SNP Replication','save',0,'outfile', ...
%             [theGroups{jj},'_',shortTitles{ii}], ...
%             'plotNumber',length(theGroups), 'plotCounter', jj)
%         hold on
%     end
% 
%     % hold off the figure
%     hold off
% 
%     % save the figure: adds _ManhattanPlot to filenames
%     name = [theGroups{jj},'_',shortTitles{ii},'_ManhattanPlot.fig'];
% 
%     % try to save the figure
%     fprintf('\nNow saving the figure \n')
%     savefig(gcf,name,'compact');

%     % remove the nan groups in both data, and the variable name to p-values
%     % and save to excel
%     bothData( isnan(bothData.AFR) | isnan(bothData.EUR), :) = [] ;


%     % make data to be used for making a scatter plot
%     plotData = bothData ;
%     plotData.logP_EUR = -log10(plotData.p_EUR) ;
%     plotData.logP_EUR = -log10(plotData.p_EUR) ;

    % save the complete data for ploting in tableau 
    bothData.Properties.VariableNames([7,8]) = {'p_EUR','p_AFR'};
    writetable(bothData(1:50,:),'SupplementaryFile5.xlsx', ...
        'Sheet',shortTitles{ii})

    % add the beta values to the data 
    curBetas = readtable(['beta_',myVars{ii},'.txt']) ;
    curBetas = curBetas(:, {'ID','BETA_META','BETA_AFR','BETA_EUR',...
        'SE_META','SE_AFR','SE_EUR'});

    % add the beta values to the snp datasets 
    bothData = innerjoin(bothData, curBetas,'Keys',{'ID'});

    % save the data excel 
    writetable(bothData,'Tableau_Reviewer_Data_heart.xlsx', ...
        'Sheet',shortTitles{ii})
end

% show the table
disp(sigSums)

% What about attempting to replicate the significant findings from one
% ancestry group in the other? For instance, are the significant findings
% in AFR at least at p<0.05 in the EUR?
clear curBetas

%% ************************** Review 2 Comment **************************

% What about doing local replication (replication across ancestry groups
% that takes into account differences in LD in the region)?

% ANSWERED 

%% ************************** Review 2 Comment **************************

% Ln 189: It is interesting to assess the differences in MAF across
% ancestries for the significant associations to see if an associated
% variant is common and detectable in one ancestry group and absent or rare
% (i.e. undetectable) in the other ancestry group as an explanation as to
% why an association in one population is not observed in the other.
% However, all of the variants shown in Supp Fig 2 are all common in both
% populations, meaning that this isn’t likely to explain why some
% associations are not captured in the other ancestry group. I do not
% understand why the authors go on to describe the variants that were most
% differentiated between the two groups as if they are more plausible
% candidates because of this differentiation. If it is thought that these
% variants being highly differentiated between ancestry groups makes them
% more likely to be important loci for pulmonary function measures compared
% to the other significantly associated variants, the authors need to
% explain that reasoning, as it is not clear.

% ANSWERED 

%% Plot the Venn diagrams of the Causal Variants 

% preallocate a table for the significant snps
snpSigTableCausal = table('Size',[1,5],'VariableTypes',...
    {'string','double','double','double','cellstr'}, ...
    'VariableNames',{'Variable','sigSnpsAFR','sigSnpsEUR', ...
    'commonSigSnps','listSigSnps'});

% --clump-kb 250 Physical distance threshold for clumping
% http://zzz.bwh.harvard.edu/plink/clump.shtml

% for the each variable
for ii = 1:length(myVars)
    
   fprintf('\n Finding the common SNPs plot for %s \n', figTitles{ii})
   
   % get the data in a loop
   curDataAFR  = readtable([myVars{ii},'_gwas_AFR.txt']) ;
   curDataEUR = readtable([myVars{ii},'_gwas_EUR.txt']) ;
   
   % load the the causal variants 
   picsAFR = readtable(['PICS2_',shortTitles{ii},'_AFR.txt']);
   picsEUR = readtable(['PICS2_',shortTitles{ii},'_EUR.txt']);
   
   % if there is not data in the african datasets
   if isempty(picsAFR)
        picsAFR = table('Size',[0,width(picsEUR)],'VariableNames',...
            picsEUR.Properties.VariableNames, ...
            'VariableTypes',{'double','cell','cell','double',...
            'double','double','double','double'});  
   end
   
%    % return only the snps that are also the linked snps
%    picsEUR = picsEUR(strcmp(picsEUR.Top_SNP, picsEUR.SNP), :);
%    picsAFR = picsAFR(strcmp(picsAFR.Top_SNP, picsAFR.SNP), :);
   
   % return only the causal variants with a PIC probablity greater than
   % 0.025
   picsAFR(picsAFR.PICS_probability < 0.1| ...
       isnan(picsAFR.PICS_probability), :) = [] ;
   picsEUR(picsEUR.PICS_probability < 0.1 | ...
       isnan(picsEUR.PICS_probability), :) = [] ;
   
   % get only the significant snps from both tables
   curDataAFR = curDataAFR(curDataAFR.P <= 5e-8,:) ;
   curDataEUR = curDataEUR(curDataEUR.P <= 5e-8,:) ; 
   
   % get the unique snps 
   [~,theUnique] = unique(curDataAFR.ID);
   
   % get only the unique snp 
   curDataAFR = curDataAFR(theUnique,:);
   curDataEUR = unique(curDataEUR);
   
   % now get only teh SnPS have have are in the causal sets
   curDataAFR = curDataAFR(ismember(curDataAFR.ID,picsAFR.SNP),:);
   curDataEUR = curDataEUR(ismember(curDataEUR.ID,picsEUR.SNP),:);
   
   % how many snps are significant in EUR and in black for the current
   % variable
   snpSigTableCausal.Variable(ii) = figTitles(ii);
   snpSigTableCausal.sigSnpsAFR(ii) = height(curDataAFR);
   snpSigTableCausal.sigSnpsEUR(ii) = height(curDataEUR) ;
   
   % get the common snps
   commonSnps = intersect(curDataAFR.ID, curDataEUR.ID) ;
   snpSigTableCausal.commonSigSnps(ii) = length(commonSnps);
   snpSigTableCausal.listSigSnps{ii} = commonSnps' ;
   
   % draw the venn diagram for the SNPs that are intersecting
   figure()
   venn([200,200],75,'EdgeColor','black')
   hold on 
   
   % add text to the figure: get numbers of intersects to add to the plots
   myNumText(1) = height(curDataAFR) - length(commonSnps) ;
   myNumText(2) = length(commonSnps);
   myNumText(3) = height(curDataEUR) - length(commonSnps) ;
  
   
   % here are text positions and the titles of the plots
   textPos = [0.25, 0.47 , 0.70] ;
   for jj = 1:length(textPos)
       text(textPos(jj),0.5,num2str(myNumText(jj)),'Units','normalized',...
           'FontWeight','bold','FontSize',16,'HorizontalAlignment','center')
   end
   
   % add the titles to the circle
   myGroups = {'AFR','EUR'};
   textPos = [0.33 , 0.66] ;
   for jj = 1:length(textPos)
       text(textPos(jj), 0.98 ,[shortTitles{ii}, ' ', myGroups{jj}], ...
           'Units','normalized', ...
           'FontWeight','bold','FontSize',14,'HorizontalAlignment','center')
   end
   hold off
   
   % save the figure
   saveas(gcf,[shortTitles{ii},'_intersect_causal.fig'],'fig')

   % sort the rows of the table
   curDataAFR = sortrows(curDataAFR,'P','ascend');
   curDataEUR = sortrows(curDataEUR,'P','ascend');
   
   % also save to the supplementary file 1
   writetable(curDataAFR,'SupplementaryFile1.xlsx','Sheet', ...
       [shortTitles{ii},' AFR'])
   writetable(curDataEUR,'SupplementaryFile1.xlsx','Sheet', ...
       [shortTitles{ii},' EUR'])
   
end

disp(snpSigTableCausal)

clear ii curDataEUR curDataAFR commonSnps ii  textPos myNumText

%% ************ AGAIN: Replot the Common Causal SNPs ************

% preallocate a table for the significant snps
snpSigTableCausal = table('Size',[1,5],'VariableTypes',...
    {'string','double','double','double','cellstr'}, ...
    'VariableNames',{'Variable','sigSnpsAFR','sigSnpsEUR', ...
    'commonSigSnps','listSigSnps'});

% for the each variable
for ii = 1:length(shortTitles)
    
   fprintf('\nPlotting common causal snps for %s \n',figTitles{ii})
   
   % get the data in a loop
   curDataAFR  = readtable(['sigSnp',shortTitles{ii},'_AFR.txt']) ;
   curDataEUR = readtable(['sigSnp',shortTitles{ii},'_EUR.txt']) ;
   
   % if there is not data in the african datasets
   if isempty(curDataAFR)
       curDataAFR = table('Size',[0,width(curDataEUR)],'VariableNames',...
            curDataEUR.Properties.VariableNames, ...
            'VariableTypes',{'cell','double','double','double'});  
   end
   
   % add the variable names
   curDataAFR.Properties.VariableNames = {'ID','CHR','POS','P'} ;
   curDataEUR.Properties.VariableNames = {'ID','CHR','POS','P'} ;
   
   % load the the causal variants 
   picsAFR = readtable(['PICS2_',shortTitles{ii},'_AFR.txt']);
   picsEUR = readtable(['PICS2_',shortTitles{ii},'_EUR.txt']);
   
   % if there is not data in the african datasets
   if isempty(picsAFR)
       picsAFR = table('Size',[0,width(picsEUR)],'VariableNames',...
           picsEUR.Properties.VariableNames, ...
           'VariableTypes',{'double','cell','cell','double',...
           'double','double','double','double'});
   end
   
%    % return only the snps that are also the linked snps
%    picsEUR = picsEUR(strcmp(picsEUR.Top_SNP, picsEUR.SNP), :);
%    picsAFR = picsAFR(strcmp(picsAFR.Top_SNP, picsAFR.SNP), :);
   
   % return only the causal variants with a PIC probablity greater than
   % 0.025
   picsAFR(picsAFR.PICS_probability < 0.1| ...
       isnan(picsAFR.PICS_probability), :) = [] ;
   picsEUR(picsEUR.PICS_probability < 0.1 | ...
       isnan(picsEUR.PICS_probability), :) = [] ;
   
   % get only the significant snps from both tables
   curDataAFR = curDataAFR(curDataAFR.P >= -log10(5e-8),:) ;
   curDataEUR = curDataEUR(curDataEUR.P >= -log10(5e-8),:) ; 
   
   % get the unique snps 
   [~,theUnique] = unique(curDataAFR.ID);
   
   % get only the unique snp 
   curDataAFR = curDataAFR(theUnique,:);
   curDataEUR = unique(curDataEUR);
   
   % now get only teh SnPS have have are in the causal sets
   curDataAFR = curDataAFR(ismember(curDataAFR.ID,picsAFR.SNP),:);
   curDataEUR = curDataEUR(ismember(curDataEUR.ID,picsEUR.SNP),:);

   % how many snps are significant in EUR and in black for the current
   % variable
   snpSigTableCausal.Variable(ii) = figTitles(ii);
   snpSigTableCausal.sigSnpsAFR(ii) = height(curDataAFR);
   snpSigTableCausal.sigSnpsEUR(ii) = height(curDataEUR) ;
   
   % get the common snps
   commonSnps = intersect(curDataAFR.ID, curDataEUR.ID) ;
   snpSigTableCausal.commonSigSnps(ii) = length(commonSnps);
   snpSigTableCausal.listSigSnps{ii} = commonSnps' ;
   
   % draw the venn diagram for the SNPs that are intersecting
   figure()
   venn([200,200],75,'EdgeColor','black')
   hold on 
   
   % add text to the figure: get numbers of intersects to add to the plots
   myNumText(1) = height(curDataAFR) - length(commonSnps) ;
   myNumText(2) = length(commonSnps);
   myNumText(3) = height(curDataEUR) - length(commonSnps) ;
  
   
   % here are text positions and the titles of the plots
   textPos = [0.25, 0.47 , 0.70] ;
   for jj = 1:length(textPos)
       text(textPos(jj),0.5,num2str(myNumText(jj)),'Units','normalized',...
           'FontWeight','bold','FontSize',16,'HorizontalAlignment','center')
   end
   
   % add the titles to the circle
   myGroups = {'AFR','EUR'};
   textPos = [0.33 , 0.66] ;
   for jj = 1:length(textPos)
       text(textPos(jj), 0.98 ,[shortTitles{ii}, ' ', myGroups{jj}], ...
           'Units','normalized', ...
           'FontWeight','bold','FontSize',14,'HorizontalAlignment','center')
   end
   hold off
   
   % save the figure
   saveas(gcf,[shortTitles{ii},'_intersect_causal2.fig'],'fig')
   
end

disp(snpSigTableCausal)

clear ii curDataEUR curDataAFR commonSnps ii  textPos myNumText

% % *************** STOP RUNNING THE CODE HERE ****************
% 
% return

%% ************************** Review 2 comments **********************

% Thanks for the nice web resource. But it seems the authors still have not
% done any serious colocalization analysis or LD checking when reporting
% signals and defining novelty. I doubt some of the reported 817
% variants/signals are actually representing the same signals, i.e.
% variants with different SNP IDs not necessary mean they are independent
% signals.

% load the LD data 
ldData = readtable('GWAScat Small PICS2 2021-09-24.txt') ;

% get the variants that are associated with the three pulmonary traits and
% return only the required rows
snpLd = readtable('snp_freq_lung.csv') ;

fprintf('\n Processing the LD data \n')
% get only the required variables
ldData = ldData(:,{'IndexSNP','LinkedSNP','x_hg19_','Rsquare'});
head(ldData) % for debugging purposes
ldData(ldData.Rsquare < 0.4, :) = [] ;

% get only the data that exist for africans and europeans
ldData = ldData( ismember(ldData.IndexSNP, snpLd.Variant) | ...
    ismember(ldData.LinkedSNP, snpLd.Variant),:) ;

% clean up the chromosome name 
ldData.Properties.VariableNames(3) = "Chrom" ;
ldData.Position = extractBefore(extractAfter(ldData.Chrom,':'),')');
ldData.Chrom = extractAfter(extractBefore(ldData.Chrom,':'), '(');

% convert to double 
ldData.Chrom = str2double(ldData.Chrom);
ldData.Position = str2double(ldData.Position); 

% remove the data with no linked snp
ldData(~contains(ldData.LinkedSNP,'rs'),:) = [] ;
head(ldData) % for debugging purposes

% save the originial copy of the snps
snpLdOG = snpLd ;

%% Now change the SNPs in the snpLd table 

% lets have the number of LD variants and get the original variable table
% to use in a loop 
therSquares = [0.8,0.6,0.4] ;
rSquareNames = {'eight','six','four'} ;
totalsnp = zeros(1,3) ;
windwSize = 1000000 ;

for jj = 1:length(therSquares)

    % get the current ld data and the curnspl data
    curLddata = ldData(ldData.Rsquare >= therSquares(jj), :) ;
    snpLd = snpLdOG;
    allLd = [] ;
    numOfLDsnps = 0 ;

    % here are the data
    for ii = 1:height(snpLd)

        % print something to the screen
        if rem(ii,50) == 0
            fprintf('\nrunning analysis for variant #%d of %d\n', ...
                ii, height(snpLd) )
        end

        % get the current snp
        curSnp = snpLd.Variant(ii) ;

        % get the snps in close proximity with the current snp -- these are
        % snps which are within 500KB of the lead snp that is predicted
        % caused
        closeSnps = snpLd( ismember( snpLd.Position , ...
            snpLd.Position(ii)-windwSize:snpLd.Position(ii)+windwSize) & ...
            snpLd.Chrom == snpLd.Chrom(ii) , : ) ;

        % remove the current snp from the table
        closeSnps(ismember(closeSnps.Variant, curSnp), :) = [] ;

        % check if the table has data
        if isempty(closeSnps)
            continue
        end

        % find if any of those snps are in LD with the lead snp
        curLd = curLddata( ismember(curLddata.IndexSNP,curSnp) , :) ;

        % find the snps are in strong LD
        curLd = curLd( ismember( curLd.LinkedSNP ,closeSnps.Variant), :) ;
        curLd = unique(curLd) ;

        % check if the table has data
        if isempty(curLd)
            continue
        end

        % print some thing the screen
        numOfLDsnps = numOfLDsnps + height(curLd) ;
        fprintf('\n The number of LD variant identified is %d\n', ...
            numOfLDsnps)

        % create the current ld struct 
        allLd = [allLd ;curLd] ;

        % get the location of those viariants
        locVariants = ismember(snpLd.Variant,curLd.LinkedSNP);

        % change those snps in the snpLd table
        snpLd.Variant(locVariants) = curSnp ;
        snpLd.Chrom(locVariants) = snpLd.Chrom(ii);
        snpLd.Position(locVariants) = snpLd.Position(ii) ;
        snpLd.Alternative(locVariants) = snpLd.Alternative(ii);
        snpLd.Reference(locVariants) = snpLd.Reference(ii);

        % let see a small table of this
        those = snpLd(locVariants,:) ;
    end

    % add the the table of table snp and add the lds results to the strct
    totalsnp(jj) = numOfLDsnps ;
    ldSnpPerR.(rSquareNames{jj}) = allLd ;

end

% display the number of LD snps
disp(numOfLDsnps)

% get the unique snps
snpLd = unique(snpLd) ;

% plot figure of the r-square values 
therSquares = cellstr( categorical(therSquares) ) ;
therSquares = categorical(strcat('>= ', therSquares)) ;

%% 
figure()
bar(therSquares,totalsnp,'DisplayName','totalsnp')
hold on 
% edit the axis and and adjust the figure
ylabel('Number of variants')
xlabel('R^2 value')
title('\bf Common Variants Based on LD-based Mapping','FontSize',16)
set(gca,'LineWidth',1.5,'Box','off','TickDir','out','FontSize',12, ...
    'FontWeight','bold')

% ****************** VECTORISED IMPLEMENTATION *******************
% add number the top of the bar graph
text(1:length(totalsnp), flip(totalsnp),  ...
    split(num2str(flip(totalsnp))) ,'vert','bottom','horiz','center'); 
% ****************** VECTORISED IMPLEMENTATION *******************
box off
hold off

fprintf('\nThe number of unique variants is %d\n', ...
    length(unique(snpLd.Variant)) )

%% save the results to excel
writetable(snpLd,'Lung Function Variants.xlsx','Sheet','AllVariants')
writetable(snpLd(snpLd.FVC_EUR == 1, :), ...
    'Lung Function Variants.xlsx','Sheet','EUR FVC Variants')
writetable(snpLd(snpLd.FEV1_EUR == 1, :), ...
    'Lung Function Variants.xlsx','Sheet','EUR FEV1 Variants')
writetable(snpLd(snpLd.PEF_EUR == 1, :), ...
    'Lung Function Variants.xlsx','Sheet','EUR PEF Variants')

clear ii curLd closeSnps

%% Compare SNPs frequency with SNPs Assocations

% Some SNPs are associated with FEV1, FVC and PEF. Are they found at the
% same or different frequencies among EUR and black

% load the UK bio bank that has SNP frequencies

% change this the snps that are present in the UK Biobank from the data
% which also included the imputed snps 
snpFreq = readtable('full_variant_qc_metrics.txt');

% return only the required table variables
snpFreq = snpFreq(:,{'chrom','pos','alt','ref','rsid','af_EUR','af_AFR'});

% change the variables names to those present in the snp comparison table
snpFreq.Properties.VariableNames = {'Chrom','Position','Alternative',...
    'Reference','Variant','EUR','AFR'};	

% add the variable to the table in a loop
for ii = 1:length(shortTitles)

    % read in the data of the current condition in
    curAFR = readtable('SupplementaryFile1.xlsx','Sheet', ...
        [shortTitles{ii},' AFR']) ;
    curEUR =  readtable('SupplementaryFile1.xlsx','Sheet', ...
        [shortTitles{ii},' EUR']) ;

    % if africa has no data 
    if isempty(curAFR)
        curAFR = table('Size',[0,width(curEUR)],'VariableNames',...
            curEUR.Properties.VariableNames, ...
            'VariableTypes',{'double','double','cell','cell',...
            'double','cell','cell'}) ;
    end

    % add those condition to the table where the SNPs is true
    locAFR = ismember(snpFreq.Variant, curAFR.ID) ;
    locEUR = ismember(snpFreq.Variant, curEUR.ID) ;

    % add the variable to the table
    snpFreq.([shortTitles{ii},'_AFR']) = locAFR ;
    snpFreq.([shortTitles{ii},'_EUR']) = locEUR ;

end

% return only the SNPs that are significant in the GWAS of the three
% variables
snpFreq = snpFreq(any(snpFreq{:,8:end},2),:); 

% add a new variable at the end of the column to show which SNPs related to
% EUR and those related to AFR
snpFreq.SigEUR = any( ...
    snpFreq{:,{'FVC_EUR','FEV1_EUR','PEF_EUR'}},2) ;

% for the africans 
snpFreq.SigAFR = any( ...
    snpFreq{:,{'FVC_AFR','FEV1_AFR','PEF_AFR'}},2) ;

% get instance where the gwas is significant in EUR
snpFreq.SigGroup(snpFreq.SigEUR==true) = {'EUR'} ;
snpFreq.SigGroup(snpFreq.SigAFR==true) = {'AFR'} ;

% get the instances where both snps are significant
bothSig = all([snpFreq.SigEUR,snpFreq.SigAFR],2) ;
snpFreq.SigGroup(bothSig)= {'Both'};


% also add for both groups 
snpFreq.SigGroup  = categorical(snpFreq.SigGroup );

% save the data to a csv file
writetable(snpFreq,'snp_freq_lung.csv')
writetable(snpFreq,'SupplementaryFile2.csv')

% disply the signficant snps for each group
summary(snpFreq.SigGroup)

% sort the rows of the data
snpFreq = sortrows(snpFreq,'EUR','descend');

clear locAFR locEUR ii curAFR curEUR 

%% Apply a chi square test to the SNP frequencies

chiTable = snpFreq(:,{'Chrom','Position','Alternative','Reference',...
    'Variant','EUR','AFR','SigGroup'}) ;

% run the data in a loop
for ii = 1:height(chiTable)

    % print something to the screen
    if rem(ii,50) == 0
        fprintf('\nRunning chi-square test # %d of %d Variables \n',...
            ii, height(chiTable))
    end
    
    % convert the variables to a categorical array and if it
    % fails then continue
    curData = chiTable(ii,{'AFR','EUR'}) ;
    curData = curData{1,:}' ;
    
    % convert the nan to zero
    if any(isnan(curData))
        curData(isnan(curData)) = 0;
    end
    
    % get the data for the current fisher exact test 
    fisherData = table('Size',[2,2],'VariableTypes',{'double','double'},...
        'VariableNames',{'snp','Nosnp'} ) ;
    
    % add the snp to the fisher exact test table and delete the count
    % variable that is not need for the fisher exact test
    fisherData.snp = round([383471;5978].*curData)  ;
    fisherData.Nosnp = round([383471;5978] - fisherData.snp) ;
    
    % now perform the chi-square test
    [~, p, stats ] = fishertest(fisherData);
    
    % add the results to the table
    % add the variable name to the table
    chiTable.pValue(ii) = p ; 
    chiTable.OddRatio(ii) = stats.OddsRatio ; 
    chiTable.LowerBound(ii) = stats.ConfidenceInterval(1);
    chiTable.UpperBound(ii) = stats.ConfidenceInterval(2);
    
end  

% add the variable names to the table and then remove the rows with missing
% results and add the adjusted p value to the table 
chiResults = addvars(chiTable, ...
    mafdr(chiTable.pValue,'BHFDR',true),'After','pValue', ...
    'NewVariableNames','FDR') ;
chiResults = sortrows(chiResults,'pValue','ascend');

% save the results to excel 
writetable(chiResults,'SNP Comparisons Imputed.csv')

clear ii curData tbl chi2 chiPvalue labels ci p stats meanAfricans  ...
    meanWhites chiTable ttestTable chiData xValues groups plotName ...
    curCol locFreq africaChiTable

%% Get the Whose Frequencies varies the most 

% get only the sig SNPs
mostSigSNPfreq = chiResults(:,{'Variant','EUR','AFR','SigGroup',...
    'pValue','FDR','OddRatio','LowerBound','UpperBound'}) ;

% get the unique snps
[~,theUnique] = unique(mostSigSNPfreq.Variant);
mostSigSNPfreq = mostSigSNPfreq(theUnique, :);

% find the difference in the snp freq
mostSigSNPfreq = addvars(mostSigSNPfreq, ...
    mostSigSNPfreq.EUR - mostSigSNPfreq.AFR,'NewVariableNames',...
    "freqDiff",'After','AFR') ;

% sort the table 
mostSigSNPfreq = sortrows(mostSigSNPfreq,'freqDiff','descend');

% plot two figures of the differences get the data and sort according the
% frequency of the snps in EUR 
plotData = mostSigSNPfreq([1:8,end-8:end],{'Variant','EUR','AFR',...
    'freqDiff','pValue'});
plotData.Variant = categorical(plotData.Variant) ;
plotData = sortrows(plotData,'freqDiff','ascend');

% produce a bar graph
figure()
bar1 = bar(plotData.Variant, plotData{:,{'EUR','AFR'}});

% change the properties of the bars graphs
set(bar1(2),'FaceColor',groupColors(1,:),'EdgeColor',groupColors(1,:));
set(bar1(1),'FaceColor',groupColors(2,:),'EdgeColor',groupColors(2,:));

set(gca,'FontSize',14,'LineWidth',2,'Box','off','TickDir','out')
ylabel('SNP Frequency')
title('SNP Frequency Comparison','FontSize',16,'FontWeight','bold')
legend({'EUR','AFR'})

% save the figure
saveas(gcf,'SNP_frequency_comparison.fig','fig')

clear plotData theUnique

%% Add the genes symbol to the SNPs 

% The gene affected may be the same even if the SNPs are different

% read in the SNP data with vep annotations camparising EUR and AFR
if ~exist('snpFreqVepLung.mat','file')
    % get the snp freq table 
    snpFreqVep = snpFreq ;
    
    % load the table has the information on the variants
    variant_qc = readtable('variant_qc_metrics.txt');
    variant_qc = variant_qc(:,{'rsid','nearest_genes'}) ;
    variant_qc.Properties.VariableNames = ["Variant","HugoSymbol"];
    
    % merge the two tables
    snpFreqVep = innerjoin(snpFreqVep,variant_qc,'Key','Variant');
    
    % move the hugosymbol to a different location 
    snpFreqVep = movevars(snpFreqVep,"HugoSymbol","After","Variant");
    
    save('snpFreqVepLung.mat','snpFreqVep');
else
    load snpFreqVepLung.mat
end

% clear some of variables
clear ans vep theUnique biomart theMissing

%% Plot the Venn Diagram of the Genes in Each SNP Set

% convert the sig group to categorical
snpFreqVep.SigGroup = categorical(snpFreqVep.SigGroup);

% get the SNPs unique to each group
EUROnlyGenes = snpFreqVep.HugoSymbol(snpFreqVep.SigGroup == "EUR");
AFROnlyGenes = snpFreqVep.HugoSymbol(snpFreqVep.SigGroup == "AFR");

% some of the snps are found near many genes so make those into a gene list
for ii = 1:length(AFROnlyGenes)
    
    % check if the genes have a comma
    if contains(AFROnlyGenes(ii),',')
        % split the multiple genes and add them to end of the cell array
        AFROnlyGenes = [AFROnlyGenes ;split(AFROnlyGenes(ii),',')] ;
    end
end
% remove the rows that have the multiple genes  and return only the
% unique genes
AFROnlyGenes(contains(AFROnlyGenes,',')) = [];
AFROnlyGenes = unique(AFROnlyGenes);

% some of the snps are found near many genes so make those into a gene list
for ii = 1:length(EUROnlyGenes)
    
    % check if the genes have a comma
    if contains(EUROnlyGenes(ii),',')
        % split the multiple genes and add them to end of the cell array
        EUROnlyGenes = [EUROnlyGenes ; split(EUROnlyGenes(ii),',')] ;
    end
end
% remove the rows that have the multiple genes  and return only the
% unique genes
EUROnlyGenes(contains(EUROnlyGenes,',')) = [];
EUROnlyGenes = unique(EUROnlyGenes);

% remove the - from all the genes 
EUROnlyGenes(ismember(EUROnlyGenes, '-')) = [] ;
AFROnlyGenes(ismember(AFROnlyGenes,'-')) = [] ;

% get only the unique genes 
EUROnlyGenes = unique(EUROnlyGenes) ;
AFROnlyGenes = unique(AFROnlyGenes) ;

% plot a venn diagram of the snps draw the venn diagram for the SNPs that
% are intersecting
figure()
venn([200,200],75,'EdgeColor','black')
hold on

% get the common genes
commonGenes = intersect(EUROnlyGenes, AFROnlyGenes ) ;

% get the likely SNPs which are those located in the same gene but they
% different 
candidateCommonGenes = snpFreqVep(ismember(snpFreqVep.HugoSymbol,...
    commonGenes),:) ;

% add text to the figure: get numbers of intersects to add to the plots
myNumText(1) = height(AFROnlyGenes) - length(commonGenes) ;
myNumText(2) = length(commonGenes);
myNumText(3) = height(EUROnlyGenes) - length(commonGenes) ;

% here are text positions and the titles of the plots
textPos = [0.25, 0.47 , 0.70] ;
for jj = 1:length(textPos)
    text(textPos(jj), 0.5 ,num2str(myNumText(jj)),'Units','normalized', ...
        'FontWeight','bold','FontSize',16,'HorizontalAlignment','center')
end

% add the titles to the circle
myGroups = {'AFR','EUR'};
textPos = [0.33 , 0.66] ;
for jj = 1:length(textPos)
    text(textPos(jj), 0.98 ,[myGroups{jj}],'Units','normalized', ...
        'FontWeight','bold','FontSize',14,'HorizontalAlignment','center')
end
hold off

% convert to to a table 
EUROnlyGenes = array2table(EUROnlyGenes,'VariableNames',...
    {'EUR SNP Genes'}) ;
AFROnlyGenes = array2table(AFROnlyGenes,'VariableNames',...
    {'AFR SNP Genes'}) ;
bothGroupsSNPs = array2table(commonGenes,'VariableNames',...
    {'Both Groups SNP Genes'} );

% save the information of the candidate common SNPs to Excel 
writetable(candidateCommonGenes,'candidateCommonGenes.xlsx') ;

% also save to the supplementray fiel 
EUROnlyGenes(ismissing(EUROnlyGenes),:) = [] ;

% save the data to excel for pathway analysis
writetable(bothGroupsSNPs,'SupplementaryFile2.xlsx',...
    'Sheet','Location of SNPs','Range','C1')
writetable(EUROnlyGenes,'SupplementaryFile2.xlsx',...
    'Sheet','Location of SNPs','Range','A1') 
writetable(AFROnlyGenes,'SupplementaryFile2.xlsx',...
    'Sheet','Location of SNPs','Range','B1')

writetable(snpFreqVep,'SupplementaryFile2.xlsx',...
    'Sheet','SNP Freq - Annon')

clear jj textPos myNumText commonGenes

%% Get the only Novel SNPs

% get the SNPs unique to each group
novelOnly =  readtable('SupplementaryFile4.xlsx','Sheet','Novel') ;
novelOnly.SigGroup = categorical(novelOnly.SigGroup);
novelOnly = unique(novelOnly(novelOnly.SigGroup == "EUR",{'HugoSymbol'}));
novelOnly = novelOnly.HugoSymbol ;

% some of the snps are found near many genes so make those into a gene list
for ii = 1:height(novelOnly)
    
    % check if the genes have a comma
    if contains(novelOnly(ii),',')
        % split the multiple genes and add them to end of the cell array
        novelOnly = [novelOnly ;split(novelOnly(ii),',')] ;
    end
end

% remove the bad genes 
novelOnly(contains(novelOnly,',')) = [] ;

writecell(novelOnly,'Novel Genes EUR.xlsx') 

%% Change the threshold p-values to a suggestive values and plot Venns

% for the each variable
for ii = 1:length(myVars)
    
    % have another loop to compare only the snps that we found across all
    % the groups significant
    for kk = 1 % for different snps
        
        fprintf('\n Finding the common SNPs plot for %s \n',figTitles{ii})
        
        % get the data in a loop
        curDataAFR = readtable([myVars{ii}, '_gwas_AFR.txt']) ;
        curDataEUR = readtable([myVars{ii}, '_gwas_EUR.txt']) ;
        
        % get only the significant snps from both tables
        if kk == 1
            curDataAFR = curDataAFR(curDataAFR.P < 1e-6,:) ;
            curDataEUR = curDataEUR(curDataEUR.P < 1e-6,:) ;
        else
            % return only the signifiacnt SNPs and correct the p-values
            % only for those SNPs 
            curDataAFR = curDataAFR(...
                ismember(curDataAFR.ID, snpFreq.Variant ),:) ;
            curDataEUR = curDataEUR( ...
                ismember(curDataEUR.ID,snpFreq.Variant),:) ;
            
            % multiply the p-values by the number of test 
            curDataAFR.P = curDataAFR.P*( ...
                length(curDataAFR.P)) ;
            curDataEUR.P = curDataEUR.P*( ...
                length(curDataEUR.P)) ;
            
            % now return only the signficant snps 
            curDataAFR = curDataAFR(curDataAFR.P < 0.05,:);
            curDataEUR = curDataEUR(curDataEUR.P < 0.05,:);
        end
        
        % get the common snps
        commonSnps = intersect(curDataAFR.ID, curDataEUR.ID) ;
        
        % draw the venn diagram for the SNPs that are intersecting
        figure()
        venn([200,200],75,'EdgeColor','black')
        hold on
        
        % add text to the figure: get numbers of intersects to add to the
        % plots
        myNumText(1) = height(curDataAFR) - length(commonSnps) ;
        myNumText(2) = length(commonSnps);
        myNumText(3) = height(curDataEUR) - length(commonSnps) ;
        
        % here are text positions and the titles of the plots
        textPos = [0.25, 0.47 , 0.70] ;
        for jj = 1:length(textPos)
            text(textPos(jj), 0.5 ,num2str(myNumText(jj)), ...
                'Units','normalized','FontWeight','bold',...
                'FontSize',16,'HorizontalAlignment','center')
        end
        
        % add the titles to the circle
        myGroups = {'AFR','EUR'};
        textPos = [0.33 , 0.66] ;
        for jj = 1:length(textPos)
            text(textPos(jj), 0.98 ,[shortTitles{ii}, ' ', ...
                myGroups{jj}],'Units','normalized', ...
                'FontWeight','bold','FontSize',14,...
                'HorizontalAlignment','center')
        end
        hold off
        
        % save the figure
        saveas(gcf,[shortTitles{ii},'_intersect_causal_1_e6.fig'],'fig')
        
        % save the snps of AFR for the lesss stringent p values
        writetable(curDataAFR,'sig_snps_suggestive.xlsx',...
            'Sheet',['AFR ', shortTitles{ii}]) ;
        writetable(curDataEUR,'sig_snps_suggestive.xlsx',...
            'Sheet',['EUR ', shortTitles{ii}]) ;
    end 
end

clear ii curDataEUR curDataAFR commonSnps ii  textPos myNumText ...
    commonSnps jj kk

%% Create a file of the location of the suggestive genes

% this will be used for enrichment analysis 
myGroups = {'AFR','EUR'};

% loop over the ethics groups   
for jj = 1:length(myGroups)
    
    % preallocate the genes 
    theGeneLoc = [];
    
    % then the variables 
    for ii = 1:length(shortTitles)
        % get the current data
        curGenes = readtable('sig_snps_suggestive.xlsx', ...
            'Sheet',[myGroups{jj}, ' ', shortTitles{ii}]) ;
        
        % add to the gene loc 
        theGeneLoc = [theGeneLoc; curGenes ];
    end
    
    % clean up the variables and return only the unique genes
    theGeneLoc = theGeneLoc.NEAREST_GENES;
    
    % some of the snps are found near many genes so make those into a gene
    % list
    for ii = 1:length(theGeneLoc)
        
        % check if the genes have a comma
        if contains(theGeneLoc(ii),',')
            % split the multiple genes and add them to end of the cell
            % array
            theGeneLoc = [theGeneLoc ; split(theGeneLoc(ii),',')] ;
        end
    end
    % remove the rows that have the multiple genes  and return only the
    % unique genes
    theGeneLoc(contains(theGeneLoc,',')) = [];
    theGeneLoc = unique(theGeneLoc);
    
    % remove the - from all the genes
    theGeneLoc(ismember(theGeneLoc, '-')) = [] ;
    
    % get only the unique genes
    theGeneLoc = unique(theGeneLoc) ;
    
    % save to excell
    writecell(theGeneLoc,'Suggistive Enrichment Genes.xlsx',...
        'Sheet',myGroups{jj} )
end

clear theGeneLoc ii jj curGenes

%% Produce a scatter plot showing the Freq Difference in SNPs

% % get the SNPs for black and the EUR
% blacksSigSNPs = snpFreq(snpFreq.SigGroup == "AFR",:) ;
% EURSigSNPs = snpFreq(snpFreq.SigGroup == "EUR",:) ;
% 
% % plot the figure
% figure()
% hold on 
% 
% % plot for the EUR
% scatter(EURSigSNPs.AFR, EURSigSNPs.EUR,10,'filled', ...
%      'MarkerFaceColor',groupColors(2,:),'Marker','o',...
%      'MarkerEdgeColor',groupColors(2,:),'MarkerFaceAlpha',0.7)
%     
% % for the AFR
% scatter(blacksSigSNPs.AFR, blacksSigSNPs.EUR,70,'filled', ...
%      'MarkerFaceColor',groupColors(1,:),'Marker','o',...
%      'MarkerEdgeColor',groupColors(1,:))
%     
% % plot the non significant data 
% nonSigFisher = snpFreq.pValue > 1e-8 ;
% scatter(snpFreq.AFR(nonSigFisher), snpFreq.EUR(nonSigFisher), ...
%     30 ,'filled','MarkerFaceColor',[ 0.5 0.5 0.5 ],'Marker','o',...
%     'MarkerEdgeColor',[ 0.5 0.5 0.5 ],'MarkerFaceAlpha',0.7)
%           
% % edit the chart elements         
% set(gca,'FontSize',14,'LineWidth',1.5,'Box','off',...
%     'XLim',[-0.03, 0.526],'YLim',[-0.03, 0.526] )
% ylabel('SNP Frequency in EUR')
% xlabel('SNP Frequency in AFR');
% legend({'GWAS EUR','GWAS AFR','No Frequency Bias'},...
%     'Location','best')
% 
% hold off

%% Identify the pathways involved

% set up the genes
try
    enrichrVars = {'GO_Molecular_Function_2021_table',...
        'GWAS_Catalog_2019_table','Reactome_2016_table',...
        'UK_Biobank_GWAS_v1_table','Elsevier_Pathway_Collection_table'} ;
    
    % set up the plot names
    plotNames = {'GO Molecular Function','GWAS Catalog','Reactome',...
        'UK Biobank GWAS','Elsevier Pathways'} ;
    
    % plot the graphs in a loop
    for ii = 1:length(enrichrVars)
        
        % % load the data
        enrichEUR = readtable([enrichrVars{ii}, ' - EUR.txt'],...
            'Format','auto');
        enrichAFR = readtable([enrichrVars{ii}, ' - AFR.txt'],...
            'Format','auto');
        
        % clean up the enriched term names
        enrichEUR.Term = strtrim( regexprep(enrichEUR.Term, ...
            {'\(+\w*','\:+\w*',')',...
            '\R-HSA-+\w*','Homo sapiens','\w* raw'},''));
        enrichAFR.Term = strtrim( regexprep(enrichAFR.Term, ...
            {'\(+\w*','\:+\w*',')',...
            '\R-HSA-+\w*','Homo sapiens', '\w* raw'},''));
        
        % create a plot and also produce the go biological data to use for
        % the venn diagrams
        enrichr_data_plotter(enrichAFR, enrichEUR,plotNames{ii} ,...
            groupColors , {'AFR','EUR'} ) ;
%         enrichrPlot_pValue(enrichAFR, enrichEUR,plotNames{ii} ,...
%             groupColors , {'AFR','EUR'} );
       
        % save to a supplementary file
        writetable(enrichEUR,'Supplementary File 3.xlsx',...
            'Sheet',[plotNames{ii},'-EUR'])
        writetable(enrichAFR,'Supplementary File 3.xlsx',...
            'Sheet',[plotNames{ii},'-AFR'])
        
    end
    
    clear enrichEUR  enrichAFR ii plotNames enrichVars
    
catch
    
end

%% Enrichments for Novel SNPs

% set up the genes
try
    
    % % load the data
    enrich1 = readtable('Novel Genes DisGeNET_table.txt',...
        'Format','auto');
    enrich2 = readtable(...
        'Novel Genes PhenGenI_Association_2021_table.txt',...
        'Format','auto');
    
    % clean up the enriched term names
    enrich1.Term = strtrim( regexprep(enrich1.Term, ...
        {'\(+\w*','\:+\w*',')',...
        '\R-HSA-+\w*','Homo sapiens','\w* raw'},''));
    enrich2.Term = strtrim( regexprep(enrich2.Term, ...
        {'\(+\w*','\:+\w*',')',...
        '\R-HSA-+\w*','Homo sapiens', '\w* raw'},''));
    
%     % create am enrichr plot
%     enrichr_data_plotter(enrich2, enrich1,'Enrichement ' ,...
%         [0.64,0.08,0.18; 0.64,0.08,0.18]  , {'DisGeNET','PhenGenI'}) ;
    
     enrichrPlot_pValue(enrich2, enrich1,'Enrichement ' ,...
     [0.64,0.08,0.18; 0.64,0.08,0.18]  , {'DisGeNET','PhenGenI'}) ;
    
    clear enrich2 enrich1
catch
end

%% Produce a scatter plot showing the Freq Difference in PICs SNPs 

% "candidateCommonGenes" will be important here 
causalSnpsEUR = [] ;
causalSnpsAFR = [] ;

% for the top 10 snps for each phenotype ;
topTop10EUR = [] ;
topTop10Blacks = [] ;

% get the SNPs in a loops 
for ii = 1:length(shortTitles)
    
  % get the current data for EUR
  causalEUR = readtable('SupplementaryFile1.xlsx', ...
      'Sheet',[shortTitles{ii},' EUR']) ;
  
  % get the current for blacks 
  causalAFR = readtable('SupplementaryFile1.xlsx',...
      'Sheet',[shortTitles{ii},' AFR']) ;
  
  % if there is not data in the african datasets
  if isempty(causalAFR)
      causalAFR = table('Size',[0,width(causalEUR)],'VariableNames',...
          causalEUR.Properties.VariableNames, ...
          'VariableTypes',{'double','double','cell','cell','double',...
          'cell','cell'});
  end
   
  % return only the causual snp from the table of snps
  causalEUR = snpFreq(ismember(snpFreq.Variant,causalEUR.ID),:) ;
  causalAFR = snpFreq(ismember(snpFreq.Variant,causalAFR.ID),:) ;

  % plot the graph before removing the causal SNPs
  figure()
  hold on
  
  % plot for the EUR
  scatter(causalEUR.AFR, causalEUR.EUR,10,'filled', ...
      'MarkerFaceColor',groupColors(2,:),'Marker','o',...
      'MarkerEdgeColor',groupColors(2,:),'MarkerFaceAlpha',0.7)
  
  % for the AFR
  scatter(causalAFR.AFR, causalAFR.EUR,70,'filled', ...
      'MarkerFaceColor',groupColors(1,:),'Marker','o',...
      'MarkerEdgeColor',groupColors(1,:))
  
%   % plot the non significant data
%   nonSigFisher = snpFreq.pValue > 1e-8 ;
%   scatter(snpFreq.AFR(nonSigFisher), snpFreq.EUR(nonSigFisher), ...
%       30 ,'filled','MarkerFaceColor',[ 0.5 0.5 0.5 ],'Marker','o',...
%       'MarkerEdgeColor',[ 0.5 0.5 0.5 ],'MarkerFaceAlpha',0.7)
  
  % edit the chart elements
  set(gca,'FontSize',14,'LineWidth',1.5,'Box','off',...
      'XLim',[-0.03, 0.526],'YLim',[-0.03, 0.526] )
  ylabel('SNP Frequency in EUR')
  xlabel('SNP Frequency in AFR');
  legend({'GWAS EUR','GWAS AFR'},'Location','best')
  title(['Causal SNPs - ',shortTitles{ii}],'FontSize',16,...
      'FontWeight','bold')
  
  hold off 
  
  % save the figure
  saveas(gca,['Causal SNPs - ', shortTitles{ii},'.png'],'png')
  
  % add the causal snps to the table 
  causalSnpsEUR = [causalSnpsEUR ; causalEUR] ;
  causalSnpsAFR = [causalSnpsAFR ; causalAFR] ;
  
end

% convert to categorical for easy accesing 
writetable(causalSnpsEUR,'Casual SNPs.xlsx','Sheet','EUR')
writetable(causalSnpsAFR,'Casual SNPs.xlsx','Sheet','AFR')

clear nonSigFisher EURSigSNPs AFRSigSNPs locLDsnp toGoLDsnps ...
    posCausal ldlength

%% *************************** Reviewer Comment **************************
% ************************************************************************

fprintf('\nPreparing file for Reviewer 3 snp comparison\n')

% Perhaps a plot comparing the p-values in each ancestry group for all
% variants that were above a certain threshold of statistical significance?
% What about attempting to replicate the significant findings from one
% ancestry group in the other?

% What about doing local replication (replication across ancestry groups
% that takes into account differences in LD in the region)?

% convert to categorical for easy accesing
causalEURComp = readtable('Casual SNPs.xlsx','Sheet','EUR');
causalAFRComp = readtable('Casual SNPs.xlsx','Sheet','AFR');

% here are the short names
shortTitles = {'FEV1','FVC','PEF'} ;

for ii = 1:length(ukGwasFiles)
    % load the gwas files in a loop 
    curAFR_file = readtable(replace(ukGwasFiles{ii},'.txt','_AFR.txt'));
    curEUR_file = readtable(replace(ukGwasFiles{ii},'.txt','_EUR.txt'));

    % rename the 5th and 6th variable which are the p values and ID
    curAFR_file.Properties.VariableNames([5,6]) = ...
        {[shortTitles{ii},'_P_AFR'],'Variant'} ;
    curEUR_file.Properties.VariableNames([5,6]) = ...
        {[shortTitles{ii},'_P_EUR'],'Variant'} ;

    % get only those column of the table 
    curAFR_file = curAFR_file(:,5:6) ;
    curEUR_file = curEUR_file(:,5:6) ;

    % add the pvalues from the AFR table 
    causalAFRComp = innerjoin(causalAFRComp,curAFR_file,'Key','Variant');
    causalEURComp = innerjoin(causalEURComp,curAFR_file,'Key','Variant');

    % also add the values from the EUR table 
    causalAFRComp = innerjoin(causalAFRComp,curEUR_file,'Key','Variant');
    causalEURComp = innerjoin(causalEURComp,curEUR_file,'Key','Variant');

end

% save the data excel
writetable(causalAFRComp,'Reviewer3Results.xlsx','AFR')
writetable(causalEURComp,'Reviewer3Results.xlsx','EUR')

%% Which Among the SNPs are eQTLs 

% This will be done after I finish with the VEP analysis

% There are difference in which genes are eQTLs depending on the tissues.
% Let's check for the this explicitly

% define the tissue types available that are aviable in the folder
if exist('/scratch/snkmus003/ukbiobank/gtexData','dir')
    addpath('/scratch/snkmus003/ukbiobank/gtexData/GTEx_Analysis_v8_eQTL')
    cd('/scratch/snkmus003/ukbiobank/gtexData/GTEx_Analysis_v8_eQTL')
    listing = ...
        dir('/scratch/snkmus003/ukbiobank/gtexData/GTEx_Analysis_v8_eQTL');
    listing = struct2table(listing) ;
    listing = listing.name(contains(listing.name,'egenes')) ;
else
    listing = dir('/Users/sinkala/Documents/MATLAB/PostDoc Analysis/GTEx_Analysis_v8_eQTL');
    listing = struct2table(listing) ;
    listing = listing.name(contains(listing.name,'egenes')) ;
end

% get the tissue names
tissues = extractBefore(listing,'.') ;

% Process the eQTL data with only p-values for SNP that are eqtls
fprintf('\n Adding eQTLs to the SNP data \n')

% run the parfor loop
for ii = 1:length(tissues)
    
    fprintf('\n Adding eQTLs of %s the SNP data \n',tissues{ii})
    
    % load the GTEx data of eQTLs
    eqtls = readtable([tissues{ii},'.v8.egenes.txt']) ;
    
    % return only the genes and the rsID and p values
    eqtls = eqtls(:,{'rs_id_dbSNP151_GRCh38p7','pval_beta'}) ;

    % EVEN HERE CONSIDER LD WHEN I ENDIFIYING EQTLS
    
    % save a copy of all the SNP that are present in the population
    eqtls.Properties.VariableNames = ['Variant',tissues(ii)] ;
    eqtls = unique(eqtls);
    
    % return only the snps that are present in the small data 
    eqtls = innerjoin(eqtls,snpFreqVep(:,{'Variant'}) );
    
    % add to the growing tables
    snpFreqVep = outerjoin(snpFreqVep,eqtls,"MergeKeys",true) ;
    
    % also get the unique snps because i need accurate statistics
    [~, theUnique ] = unique(snpFreqVep.Variant,'first') ;
    snpFreqVep = snpFreqVep(theUnique,:) ;
end

% remove the nan data from 
snpFreqVep(isnan(snpFreqVep.Position), :) = [] ;
snpFreqVep = sortrows(snpFreqVep,'Adipose_Subcutaneous','ascend');

try
    cd('/scratch/snkmus003/ukbiobank/lungAnalysis')
catch
end

%% Clustergram of the eQTLs 

% The eQTLs should be transposed such that the genes are on top, followed
% by the ethic groups.

% get the SNPs are available in Gtex 
snpGtex = snpFreqVep ;

% remove the snps that are not eqlts 
locAdipose = find( ismember(snpGtex.Properties.VariableNames, ...
    'Adipose_Subcutaneous'), true) ;

% remove the brian tissues and other tissues that not needed 
snpGtex(:,contains(snpGtex.Properties.VariableNames, 'Brain_')) = [] ;
snpGtex = removevars(snpGtex ,{'Vagina','Uterus','Testis',...
    'Skin_Sun_Exposed_Lower_leg','Skin_Not_Sun_Exposed_Suprapubic', ...
    'Ovary','Nerve_Tibial','Minor_Salivary_Gland', ...
    'Cells_EBV-transformed_lymphocytes','Breast_Mammary_Tissue',...
    'Cells_Cultured_fibroblasts','Prostate',...
    'Adipose_Subcutaneous','Pituitary'} ) ;

% convert the number to negative logarithm and those points with missing
% data to 0
snpGtex{:,locAdipose:end} = -log10(snpGtex{:,locAdipose:end}) ;

% remove the nan values 
size(snpGtex)
snpGtex( all(isnan(snpGtex{:,locAdipose:end}),2), :) = [] ;
size(snpGtex)

% make the values in the table either one or zero
eqtlsVars = snpGtex{:,locAdipose:end} ;
eqtlsVars(isnan(eqtlsVars)) =  0 ;
eqtlsVars(eqtlsVars ~=0) =  1 ;
snpGtex{:,locAdipose:end} = eqtlsVars ;

clear eqtlVars

% create labels for correalation matrix for names in which the underscore
% has been replaced with a hyphen
theVariants = snpGtex.Variant ;
clustData = snpGtex{:,locAdipose:end}' ;

% throw in an assertion 
assert(~any(all(clustData == 0)) )

% get the tissue names 
theTissues = replace( snpGtex.Properties.VariableNames(locAdipose:end),...
    '_','-');

% produce the clustergram 
cgo = clustergram(clustData,'Colormap',redbluecmap,...
    'RowLabels',theTissues,'ColumnLabels',theVariants,'Linkage',...
    'complete','ColumnPDist','correlation','Standardize','row');

% Produce a clustered heatmap of the highlight table in the previous section

% First cluster the data to get the groups that are clustered that will be
% in the heatmap 

% these drug names are as given below. Now arrage the drugs in
% the table according to this order
[~,locX] = ismember(flip(cgo.ColumnLabels), theVariants) ;
[~,locY] = ismember(cgo.RowLabels,theTissues) ;
heatData = clustData(locY,locX);

% create a multiple plot figure in MATLAB
figure(); clf
set(gcf,'position',[100,50,500,600]);

% the first number is how far the figure will be from the x-axis and the
% seceond number is now far the figure will be from the y-axis. The third
% number is was far the figure will run across the figure bar and the last
% number is far it will displaced allow the y-axis

axes('position',[0.2, 0.13, 0.60, 0.42]);
heatmap(cgo.ColumnLabels, flip(cgo.RowLabels), heatData,...
    'Colormap',summer,'ColorbarVisible','off');

% create a heatmap
[~,locX] = ismember(cgo.ColumnLabels,snpGtex.Variant) ;
snpGtex = snpGtex(locX,:);

% throw in an assertion 
assert( all(strcmp(cgo.ColumnLabels',snpGtex.Variant)) )

% creat the horizontal bar graph to put on top of the heatmap map showing
% create the heatmap this should be taken from the tcga processed data or
% the second clustergram
% Initialise the variables
toAddFirst = {'SigGroup','Chrom'} ;

barNames1 = {'Groups Affected','Chromosome'};

% define the colors 
% create colors for al the groups that are to be used in the top bar charts
classesColors.SigGroup = groupColors ;
classesColors.Chrom = rand(numel(unique(snpGtex.Chrom)),3) ;

% initialise the values for the for loopp
xInitial = 0.2 ; yInitial = 0.56 ; xEndPos = 0.60 ; ySize = 0.03 ;
increaseby = 0.04;

% now make the plot
for ii = 1:length(toAddFirst)
    % add the groups
    axes('position',[xInitial,yInitial, xEndPos, ySize]);
    hold on
    barData2 = snpGtex.(toAddFirst{ii})' ;
    if iscategorical(barData2) 
        % lengendNames =
        barData2 = double(barData2) ;
    elseif iscell(barData2)
        barData2 = str2double(barData2) ;
    end
    
    % plot the heatmap for continous data and barplot for categorical data
    ultraBars(barData2,classesColors.(toAddFirst{ii}), ...
        barNames1{ii})
        
    % add the name of the cluster to the left of bar
    dim = [0.09 yInitial 0.11 ySize];
    annotation('textbox',dim,'String',barNames1{ii},'FitBoxToText','on',...
        'FontSize',10,'FontWeight','bold','EdgeColor','none',...
        'HorizontalAlignment','right');
    
    % add the cluster names to the right of the colorbar
    hold on
    if length(unique(barData2)) < 10
        xdistance = 0.03 ;
        boxPos = [xEndPos+0.21 yInitial 0.02 0.025] ;
        boxColors = classesColors.(toAddFirst{ii}) ;
        boxNumbers = unique(barData2) ;
        for jj = 1:length(unique(barData2))
            % add a box for each values 
            annotation('rectangle',boxPos,'FaceColor',boxColors(jj,:), ...
                'EdgeColor',[1 1 1])
            annotation('textbox',boxPos,'String',...
                num2str(boxNumbers(jj)),'FitBoxToText','on',...
                'FontSize',13,'FontWeight','bold','EdgeColor','none',...
                'HorizontalAlignment','center',...
                'VerticalAlignment','middle','Color',[1 1 1]) ;
            boxPos(1) = boxPos(1) + xdistance ;
        end
    end
    
    % increase the value to change the colors and plot positions
    yInitial = yInitial + increaseby;
end

% create a dendogram for the heatmap that i will produce
axes('position',[xInitial,yInitial, xEndPos, 0.12]);
tree = linkage(clustData','average');
leafOrder = optimalleaforder(tree,pdist(clustData'));
H = dendrogram(tree,0,'Reorder',leafOrder,'ColorThreshold',10);
set(H,'LineWidth',0.5)
set(gca,'YTickLabel',[],'Visible','off');

hold off 

%% Plot comparing the frequency of variant 

% get the data of EUR and AFR 
EURData = snpFreqVep.EUR(snpFreqVep.SigGroup == "EUR") ;
AFRData = snpFreqVep.AFR(snpFreqVep.SigGroup == "AFR") ;
curData = [AFRData ; EURData ] ; 

% create a categorical array for plotting
groupsHere = [ true(length(AFRData),1); false(length(EURData),1) ];
groupsHere = renamecats( categorical(double(groupsHere)) ,{'1','0'}, ...
    {'AFR','EUR'} );

% perform a ttest with unequal variable assumed
[~,p] = ttest2(EURData,AFRData,'Vartype','unequal') ;

% plot the data
colourBoxPlot(curData,groupsHere,groupColors,true)
hold on

% annotate the plots
ylabel('SNPs Frequency')
title('\bf Distribution of GWAS Significant SNPs','FontSize',16)

% put the p value on the plot depending on its value
if p < 0.0001
    text( 0.5, 0.9 , ['p ', convertPValue2SuperScript(p)],...
        'Units','normalized','FontWeight','bold','FontSize',14, ...
        'HorizontalAlignment','center')
else
    text( 0.5, 0.9 ,sprintf('p = %0.4f',p),...
        'Units','normalized','FontWeight','bold','FontSize',14, ...
        'HorizontalAlignment','center')
end

hold off

clear curData groupsHere p ci stats 

%% What Phenotypes are Associated with the Lead PICs SNPs

% snps the casaul snps and add gwas catalogue anotations to that table 
causalSNPs = [causalSnpsEUR ; causalSnpsAFR ] ;

% get only the unique SNPs to increase computation
[~,theUnique] = unique(causalSNPs.Variant);
causalSNPs = causalSNPs(theUnique,:);

% ******************* NEW SECTION OF CODE for LD ***********************

% Load the lds scores of the snps from the GWAS catalogue based on the PIC2
% calculations and return only the rows with R-square values > 0.05
fprintf('\n Loading the PIC2 LD data for the GWAS catalog \n') 
ldData = readtable("GWAScat Small PICS2 2021-09-24.txt");

% let see when I get an error here on the cluster
head(ldData)

fprintf('\n Processing the LD data \n')
% get only the required variables
ldData = ldData(:,{'IndexSNP','LinkedSNP','x_hg19_','Rsquare'});
head(ldData) % for debugging purposes
ldData(ldData.Rsquare < 0.05, :) = [] ;

% get only the data that exist for africans and europeans
ldData = ldData(ismember(ldData.IndexSNP, causalSNPs.Variant),:) ;

% clean up the chromosome name 
ldData.Properties.VariableNames(3) = "Chrom" ;
ldData.Chrom = extractAfter(extractBefore(ldData.Chrom,':'), '(');
% ldData.Chrom = str2double(ldData.Chrom);

% remove the data with no linked snp
ldData(~contains(ldData.LinkedSNP,'rs'),:) = [] ;

% repeat the matrix and add it to the causal SNPs table
ldEUR = ldData(ismember(ldData.IndexSNP,causalSnpsEUR.Variant),...
    {'IndexSNP','LinkedSNP','Chrom','Rsquare'}) ;
ldAFR = ldData(ismember(ldData.IndexSNP,causalSnpsAFR.Variant),...
    {'IndexSNP','LinkedSNP','Chrom','Rsquare'}) ;

% Keep copies of the original data 
ldEURog = ldEUR;
ldAFRog = ldAFR ;

% reshape the data so that it can be vercat with the reported snps
ldEUR = array2table(...
    [ldEUR{:,{'IndexSNP','Chrom'}} ; ldEUR{:,{'LinkedSNP','Chrom'}}], ...
    'VariableNames',{'Variant','Chrom'}) ;
ldAFR = array2table(...
    [ldAFR{:,{'IndexSNP','Chrom'}} ; ldAFR{:,{'LinkedSNP','Chrom'}}], ...
    'VariableNames',{'Variant','Chrom'}) ;

% and a column to the last column the sig group
ldEUR = addvars(ldEUR, repmat({'EUR'},height(ldEUR),1),...
    'NewVariableNames',{'SigGroup'}) ;
ldAFR = addvars(ldAFR, repmat({'AFR'},height(ldAFR),1),...
    'NewVariableNames',{'SigGroup'}) ;

% convert to the write variable form 
ldAFR.SigGroup = categorical(ldAFR.SigGroup);
ldEUR.SigGroup = categorical(ldEUR.SigGroup);
ldEUR.Chrom = str2double(ldEUR.Chrom) ;
ldAFR.Chrom = str2double(ldAFR.Chrom) ;

% reshape the data so that it can be vercat with the reported snps
ldEURog = array2table(ldEURog{:,{'IndexSNP','LinkedSNP'}}, ...
    'VariableNames',{'IndexSNP','Variant'}) ;
ldAFRog = array2table(ldAFRog{:,{'IndexSNP','LinkedSNP'}}, ...
    'VariableNames',{'IndexSNP','Variant'}) ;

% put the variant table together 
ldOgBoth = unique([ldAFRog ;ldEURog ]) ;

% remove the snps that have the same names in both columns 
ldOgBoth(strcmp(ldOgBoth.IndexSNP,ldOgBoth.Variant), :) = [] ;

% get the unique causal snps
causalSNPs = unique(causalSNPs);

%% ********************************************************************

% use the lead and casual snps data 
% loading GWAS catalog data 
fprintf('\n Loading the GWAS catalog data \n') 
gwas = readtable('gwasCatalogue.tsv','FileType','text');
locSNPs = find(ismember(gwas.Properties.VariableNames,'SNPS'),true) ;
gwas.Properties.VariableNames(locSNPs) = "Variant" ;

% snps the casaul snps and add gwas catalogue anotations to that table 
snpsReport = causalSNPs(:,{'Variant','Chrom','SigGroup'});
snpsReport.SigGroup = categorical(snpsReport.SigGroup);

% also create a table for the causal snp freq 
causal_snpFreq = causalSNPs;

% set only the unique snps
% snpsReport = unique([snpsReport; ldEUR ;ldAFR]);

% get the previous report from teh GWAS catalogue
snpsReport = outerjoin(snpsReport,gwas,'Key','Variant','MergeKeys',true) ;

% % clean up the data 
% snpsReport(isnan(snpsReport.Chrom),:) = [] ;

% get the number of novel SNPs and the reported Snps
novelSnps = snpsReport(isnan(snpsReport.PUBMEDID),1:3);
reportedSnps = snpsReport(~isnan(snpsReport.PUBMEDID),:);
lungAssocSnps = reportedSnps(contains(reportedSnps.DISEASE_TRAIT, ...
    {'lung','Pulmonary','airway','airflow','smok','BMI','asthma',...
    'obstru','length','Body mass index','Body size','Emphysema','FEV1',...
    'Height','Hip circumference','FVC','weight','expiratory',...
    'Pneumonia','broncho','Respiratory','Waist'}, 'IgnoreCase',true), :) ;

% get the actual snps related to lung function 
lungReportedSnps = lungAssocSnps( ...
    contains(lungAssocSnps.DISEASE_TRAIT,...
    {'lung','Pulmonary','air','smok','asthma','obstru', ...
    'Emphysema','FEV1','FVC','expiratory',...
    'Pneumonia','broncho','Respiratory'}, 'IgnoreCase',true) , :) ;

% add the snps related to the waist back to the lung associatd snps
lungAssocSnps = [lungAssocSnps; lungReportedSnps( ...
    contains(lungReportedSnps.DISEASE_TRAIT,{'waist'}, ...
    'IgnoreCase',true),:)] ;
lungReportedSnps(contains(lungReportedSnps.DISEASE_TRAIT,...
    {'Lung Cancer','Adenocarcinoma','cancer','waist'},...
    'IgnoreCase',true) , :) = [] ; 

% merge with the snp freq table 
novelSnps.Properties.VariableNames(1) = "Variant" ; 
novelSnps = innerjoin(novelSnps(:,1),snpFreqVep,'Key',"Variant") ;

% return only the lung related SNPs 
lungAssocSnps(ismember( ...
    lungAssocSnps.Variant,lungReportedSnps.Variant),:) = [] ;

% load the pharos disease relations to be use to be obtain gene linked to
% disease
pharosDx = readtable('Pharos Disease Linkage.csv') ;
pharosDx = pharosDx(:,{'Symbol','DiseaseDataSource',...
    'AssociatedDiseaseP_value','AssociatedDiseaseSourceID',...
    'LinkedDisease'}) ;
pharosDx.Properties.VariableNames(1) = "HugoSymbol" ;

% get the gene related SNPs 
lungDxRelatedSnps = novelSnps(:,{'Variant','Chrom','Position',...
    'HugoSymbol','SigGroup'}) ; 

% add the disease phenotype to the data and get only disease associated
% with lung function
lungDxRelatedSnps = innerjoin(lungDxRelatedSnps,pharosDx,...
    'Key','HugoSymbol') ;
lungDxRelatedSnps = lungDxRelatedSnps( contains(...
    lungDxRelatedSnps.LinkedDisease, {'lung','Asthma','pulmonary',...
    'Airflow','bronchi','Hamman-Rich','asphyxia'},'IgnoreCase',true),:) ;

% move the lung cancer from the lungReportedSnps to the lungDxRelatedSnps 
lungCancer = lungReportedSnps( ...
    contains( lungReportedSnps.DISEASE_TRAIT, ...
    {'Lung Cancer','Adenocarcinoma','cancer'},'IgnoreCase',true), :) ;
lungCancer = innerjoin(lungCancer(:,{'Variant'}),snpFreqVep,...
    'Key',"Variant") ;
lungCancer = lungCancer(:,{'Variant','Chrom','Position',...
    'HugoSymbol','SigGroup'}) ; 
lungCancer = innerjoin(lungCancer,pharosDx,'Key','HugoSymbol') ;
lungDxRelatedSnps = [lungDxRelatedSnps ; lungCancer ] ;

% get the eQTL related SNPs that are also signficant
eQTLrelatedSnps = novelSnps(~isnan(novelSnps.Lung),...
    {'Variant','Chrom','Position','HugoSymbol','Lung','SigGroup'}) ; 
eQTLrelatedSnps(eQTLrelatedSnps.Lung > 0.05, :) = [] ;
eQTLrelatedSnps.Properties.VariableNames(5) = "Lung eQTLs pValue" ;

% remove the eQTLsnps and lungDxsnps from the novel snps table 
novelSnps(ismember(novelSnps.Variant,lungDxRelatedSnps.Variant),:) = [] ;
novelSnps(ismember(novelSnps.Variant,eQTLrelatedSnps.Variant),:) = [] ;

% get some summary statistical of the novel SNPs, related SNPs, reported
% SNPs, and SNPs related by gene function (gene in which the SNP is located
% are assoicted with lung function
snpTypeTable = table('Size',[2,6],'VariableTypes',...
    {'String','double','double','double','double','double'},...
    'VariableNames',{'Race','PulmonaryReported','PulmonaryAssociated',...
    'LungDiseaseAssociated','eQTL','Novel'}) ;
snpTypeTable.Race = {'AFR';'EUR'} ;

% get the classess of SNPs in a table also add the classes to top snp table

% get the reported snps in a loop
for ii = 1:5
    % print something to the screen 
    fprintf('\nProcessing the known and novel SNPs number: %d\n',ii)
    
    % get the current data for the loop 
    switch ii
        case 1
            curData = lungReportedSnps;
        case 2
            curData = lungAssocSnps;
        case 3
            curData = lungDxRelatedSnps;
        case 4
            curData = eQTLrelatedSnps;     
        otherwise
            curData = novelSnps;
    end
   
    % add the names to snp table and theTop4 snp table and the snpFreq
    % table
    causal_snpFreq.Evidence( ...
        ismember(causal_snpFreq.Variant,curData.Variant)) = ...
        snpTypeTable.Properties.VariableNames(ii+1) ;
    
    curData(~ismember(curData.Variant,causalSNPs.Variant), :) = [] ;
    
    % remove the lung reported snps
    if ii > 1
        % find which among those snps are member of the lungReportedSnps
        locThem = ismember(curData.Variant,lungReportedSnps.Variant);
        
        % now this is the clean set of SNPs that that is in LD with the the
        % reported snps
        curData(locThem,:) = [] ;
    end
       
    % the SNPs are repeated in the data therefore only get the unique for
    % the statisics
    [~,locUnique] = unique(curData.Variant) ;
    statData = curData(locUnique,:) ;
    statData = statData.SigGroup ;
    
    % get the number of AFR and EUR snps
    snpTypeTable.(ii+1) = countcats(statData) ;
    
    % save to excel 
    writetable(curData,'Supplementary File 4.xlsx', ...
        'Sheet',snpTypeTable.Properties.VariableNames{ii+1})
  
end

% save the top five causal snps to excel 
writetable(snpTypeTable,'theTop Causal SNPs each.xlsx','Sheet',...
    'snp Type Table');

% make the rows with missing data into Novel 
causal_snpFreq.Evidence( ...
    cellfun(@isempty,causal_snpFreq.Evidence)) = {'Novel'} ;
summary(categorical(causal_snpFreq.Evidence))

clear locSNPs ii locVep reportNames casualSNPs statData lungCancer ...
    snpsToKeep theKnownSnps ldAFRog ldEURog ldCurData toRemoveVariants

%% Create a table of the most Sig SNPs

% preallocate the top 5 snps in each groups
theTop5 = [] ;
topNumber = 10 ;

% get the top10 snps for each variable 
top10Table = causal_snpFreq(:,{'Variant','Evidence'});

% loop over the vars
for ii = 1:length(myVars)
    
    % loop over the groups
    for jj = 1:length(myGroups)
        
        fprintf('\nGetting SNPs for LD plotting data for %s\n',...
            figTitles{ii})
        % get the data in a loop
        topSNPs  = readtable('SupplementaryFile1.xlsx','Sheet', ...
            [shortTitles{ii},' ',myGroups{jj}]) ;
        
        % remove the unrequired variables
        topSNPs = removevars(topSNPs,{'REF','ALT'});
        
        % change the variable names so that I can do an inner join
        topSNPs.Properties.VariableNames = ...
            {'Chrom','Position','P_gwas','Variant','NEAREST_GENES'} ;
        
        % merge the two table % If the table has data
        if ~isempty(topSNPs)
            topSNPs = innerjoin(topSNPs,top10Table);
        else
            continue
        end
        
        % sort the data based on the chromosome
        topSNPs = sortrows(topSNPs,'Chrom','ascend');
        
        % get the top 5 SNPs in each table and add a new variable name for
        % the group and the variable
        if height(topSNPs) > topNumber
            % add the AFR data to the top 5 table if we have more than 5
            % significant SNPs for that comparison
            [~,locMostSig] = mink(topSNPs.P_gwas,topNumber) ;
            theTop5 = [ theTop5; addvars( topSNPs(locMostSig,:), ...
                repmat(myGroups(jj),length(locMostSig),1) , ...
                repmat(shortTitles(ii), length(locMostSig),1), ...
                'NewVariableNames',{'Ethnicity','Measure'} ) ] ;
        else
            % add all the snps to the table
            theTop5 = [ theTop5; addvars(topSNPs, ...
                repmat(myGroups(jj),height(topSNPs),1) , ...
                repmat(shortTitles(ii), height(topSNPs),1), ...
                'NewVariableNames',{'Ethnicity','Measure'} ) ] ;
        end
        
    end
    
    
end

% move around the variable names 
theTop5.ChromPos = cellstr([ num2str(theTop5.Chrom), ...
    repmat(':',height(theTop5),1),num2str(theTop5.Position)]);
theTop5 = theTop5(:,{'Variant','NEAREST_GENES','Ethnicity','Measure',...
    'ChromPos','P_gwas','Evidence'});

% sort the rows 
theTop5 = sortrows(theTop5,'Ethnicity','ascend') ;
theTop5 = sortrows(theTop5,'Measure','ascend') ;

% save the data to 
writetable(theTop5,'theTopSNPsTable.xlsx') ;

clear ii jj locMostSig topNumber

%% Look at the ethic groups in the gwas catalogue studie of pulmonary fxn

% get the actual snps related to lung function 
gwasLungBias = gwas( ...
    contains(gwas.DISEASE_TRAIT,...
    {'lung','Pulmonary','airway','smok','asthma','obstru', ...
    'Emphysema','FEV1','FVC','expiratory',...
    'Pneumonia','broncho','Respiratory'}, 'IgnoreCase',true) , :) ;

% add the snps related to the waist back to the lung associatd snps
gwasLungBias(contains(gwasLungBias.DISEASE_TRAIT,...
    {'Lung Cancer','Adenocarcinoma','cancer','waist','hair'},...
    'IgnoreCase',true) , :) = [] ; 

% define the EUR, AFR and others
EURtudies = {'europ','EUR','British','Finnish','French', ...
    'Dutch','Swedish','Icelandic'} ;
africaStudies = {'africa','black'} ;
othersStudies = {'Hispanic','Caribbean','Latino','Asian','Korean',...
    'Chilean','Japanese','Chinese','Peruvian','Puerto Rican',...
    'Hawaiian','Costa Rican','Bangladeshi','Indian','Hutterite','Mexican'};

% let get information of the study 
gwasLungBias = addvars( gwasLungBias,  ...
    contains(gwasLungBias.INITIALSAMPLESIZE,EURtudies,...
    'IgnoreCase',true) | ...
    contains(gwasLungBias.REPLICATIONSAMPLESIZE,EURtudies,...
    'IgnoreCase',true),...
    contains(gwasLungBias.INITIALSAMPLESIZE,africaStudies, ...
    'IgnoreCase',true) | ...
    contains(gwasLungBias.REPLICATIONSAMPLESIZE,africaStudies,...
    'IgnoreCase',true),...
    contains(gwasLungBias.INITIALSAMPLESIZE,othersStudies, ...
    'IgnoreCase',true) | ...
    contains(gwasLungBias.REPLICATIONSAMPLESIZE,othersStudies,...
    'IgnoreCase',true),...
    'After','PUBMEDID',...
    'NewVariableNames',{'Europeans','AFR','Others'} ) ;

% include Hispanic British Caribbean Latino Asian Finnish Korean  ...
% Chilean Japanese Chinese French Peruvian Puerto Rican Dutch Swedish ...
% Hawaiian Costa Rican Bangladeshi Indian Hutterite Mexican Icelandic

% get the unique link to the studies
[~, theUnique] = unique(gwasLungBias.LINK);
gwasLungBias = gwasLungBias(theUnique,:) ;
gwasLungBias = movevars( gwasLungBias,...
    {'Variant','MAPPED_GENE','REGION','CHR_ID','CHR_POS'},'Before',1) ;

% *********** do the same for the study bias for reported studies *********
studyBias = lungReportedSnps ; % lungReportedSnps

% let get information of the study 
studyBias = addvars(studyBias,  ...
    contains(studyBias.INITIALSAMPLESIZE,EURtudies,...
    'IgnoreCase',true) | ...
    contains(studyBias.REPLICATIONSAMPLESIZE,EURtudies,...
    'IgnoreCase',true),...
    contains(studyBias.INITIALSAMPLESIZE,africaStudies, ...
    'IgnoreCase',true) | ...
    contains(studyBias.REPLICATIONSAMPLESIZE,africaStudies,...
    'IgnoreCase',true),...
    contains(studyBias.INITIALSAMPLESIZE,othersStudies, ...
    'IgnoreCase',true) | ...
    contains(studyBias.REPLICATIONSAMPLESIZE,othersStudies,...
    'IgnoreCase',true),...
    'After','PUBMEDID',...
    'NewVariableNames',{'EUR','AFR','Others'} ) ;

% get the unique link to the studies
[~, theUnique] = unique(studyBias.LINK);
studyBias = studyBias(theUnique,:) ;

% remove the variables that are not requried 
studyBias = removevars(studyBias(:,1:21), {'SigGroup'}) ;

% get the number of samples for each study
gwasLungBias = addvars(gwasLungBias, gwasLungBias.INITIALSAMPLESIZE, ...
    'NewVariableNames',"SampleSize", 'After',"CHR_POS") ;

% delete the commas between the numbers 
gwasLungBias.SampleSize = replace( gwasLungBias.SampleSize,',','') ;

% now return only the numbers from the text 
for ii = 1:height(gwasLungBias)
    
    % get the location of the numbers
    sampSizeString = replace( gwasLungBias.SampleSize(ii) ,',','') ;
    numLoc = cell2mat( regexp(sampSizeString,'\d') ) ;
    
    % some times there are more than one samples 
    jumpLocs = [] ;
    
    % check for a jump in the numbers
    for jj = 2:length(numLoc)
        if numLoc(jj) - numLoc(jj-1) > 1           
            jumpLocs = [jumpLocs; jj-1 ] ;
        end   
    end
    
    % replace with the number of samples
    if isempty(jumpLocs)
        % when only one number exists
        sampSizeString = char(sampSizeString) ;
        gwasLungBias.SampleSize(ii) = cellstr(sampSizeString(numLoc)) ;
    else % when more than 1 number exists
        totalSamples = 0 ;
        
        % convert the to a char 
        sampSizeString = char(sampSizeString) ;
        
        % add the number of samples to the total sample
        if length(jumpLocs) == 1
            totalSamples = totalSamples + ...
                str2double(sampSizeString(numLoc(1:jumpLocs(1)))) ;
            totalSamples = totalSamples + ...
                str2double(sampSizeString(numLoc(jumpLocs(1)+1:end))) ;
        else          
            % add the number of samples to the total sample
            totalSamples = totalSamples + ...
                str2double(sampSizeString(numLoc(1:jumpLocs(1)))) ;
            
            % add the subsequent number to the total
            for kk = 1:length(jumpLocs)-1
                % the number of smaples
                curSamples = sampSizeString( numLoc( ...
                    jumpLocs(kk)+1:jumpLocs(kk+1) )) ;
                
                % add the number to the total sample size
                totalSamples = totalSamples + str2double(curSamples) ;
            end
            
            % add the last number
            totalSamples = totalSamples + str2double( ...
                sampSizeString(numLoc(jumpLocs(end)+1:end)) ) ;
           
        end
        % add the number of samples to the table
        gwasLungBias.SampleSize(ii) = cellstr(num2str(totalSamples)) ;
    end
end

% convert the numbers to double 
gwasLungBias.SampleSize = str2double(gwasLungBias.SampleSize);

% save to excel 
writetable(gwasLungBias,'studyBias Tableau.xlsx','Sheet','lung GWAS Bias') 
writetable(studyBias,'studyBias Tableau.xlsx','Sheet','sig GWAS Bias')

clear theUnique

%% ************************** REVIEWER COMMENT ***************************
% The increasing gap over the past 5 years would seem to coincide with
% broad accessibility of the UKB data – it this one of the drivers of this
% gap (as most analyses tend to just exclude non-European ancestry
% individuals)?

figure()
histogram(log10(gwasLungBias.SampleSize))

% adjust the figure
set(gca,'FontSize',12,'LineWidth',1,'FontSize',12 ,'Box','off',...
    'FontWeight','bold')
ylabel('Number of studies')
xlabel('Log 10 (sample size)')
title('Sample Size Across Studies','FontSize',16)

biobank = gwasLungBias( contains(gwasLungBias.STUDY, ...
    {'UK Biobank','UK','Biobank'}), :)

%% save a file of the data 

lungBias = gwasLungBias ;
lungBias = removevars(lungBias, {'REPLICATIONSAMPLESIZE',...
    'REPORTEDGENE_S_','UPSTREAM_GENE_ID','DOWNSTREAM_GENE_ID',...
    'SNP_GENE_IDS','UPSTREAM_GENE_DISTANCE','DOWNSTREAM_GENE_DISTANCE',...
    'MERGED','SNP_ID_CURRENT','INTERGENIC','PLATFORM_SNPSPASSINGQC_',...
    'CNV','PVALUE_MLOG','P_VALUE_TEXT_'});
lungBias.Properties.VariableNames{6} = 'TotalSampleSize';
lungBias = movevars(lungBias,'INITIALSAMPLESIZE','After','TotalSampleSize');

% save to excel 
writetable(lungBias,'Supplementary File 5.xlsx','Sheet',...
    'Lung GWAS Studies')

%% Let plot some PCA data 

% This file contains the information from the PCA analysis on 101284 SNPs.
% It is provided so that researchers can project their own samples onto the
% sample principal components as used within UKBiobank. The file is in MAP
% format for use with the SHELLFISH program
% (http://www.stats.ox.ac.uk/~davison/software/shellfish/shellfish.php) for
% performing PCA. Columns are tab-separated and defined as follows:

% Chromosome number; rs ID; cM location (can be anything, but the column
% must be present if the Shellfish program is used); Physical location;
% Allele 1; Allele 2; allele frequency used when computing the PCA PC1
% loading; PC2 loading; . . . PC15 loading.

pcaLoads = readtable('snp_pca_map.txt') ;

%% Do the FVC plot against BMI percetile

% A combined one and then one separate for AFR and EUR
bmiData = data(:,[{'EthnicGroupBW','BodyMassIndex_BMI_', ...
    'AgeAtRecruitment','StandingHeight'},myVars]) ;

% make the variable names easier to work with 
bmiData.Properties.VariableNames = {'Race','BMI','Age','Height',...
    'FVC','FEV1','PEF'};

% remove the Nan values from teh data 
bmiData(any(isnan(bmiData{:,2:end}),2),:) = [] ;

% add the FEV1/FVC Percentage to the table: COPD is defined as a reduced
% ratio of post-bronchodilator forced expiratory volume in 1 s (FEV1) to
% forced vital capacity (FVC) (post-bronchodilator FEV1/FVC < 0.70)
bmiData = addvars( bmiData, bmiData.FEV1 ./ bmiData.FVC , 'After', ...
    'PEF','NewVariableNames','FEV1/FVC') ;

% set up the counfounder names and make a new short title and weather to
% save the figure
myConfounders = {'BMI','Age','Height'} ;
bmiShortTitles = [shortTitles, 'FEV1/FVC'] ;
bmiVarUnits = [varUnits, '%'];
bmiFigTitles = [ shortTitles , 'FEV1/FVC'] ;
toSave = false ;

for kk = 1:length(myConfounders)
    
    % get the current confounder and the confounder variable in the table 
    curConf = myConfounders{kk} ;
    confVar = [lower(curConf),'Perc'] ;
    
    % Find the 40th and 60th percentiles of the elements of X.
    confounderPerc = prctile(bmiData.(curConf),[5:5:95]) ;
    
    % add the tenth percentile
    locPercentile = bmiData.(curConf) <= confounderPerc(1) ;
    bmiData.(confVar)(locPercentile) = 5;
    
    % preallocate the percentile to start with
    percVar = 10 ;
    
    % now add a new variable to the table called BMI percentile
    for jj = 1:length(confounderPerc)
        
        % get the location of the percentiles
        locPercentile = bmiData.(curConf) >= confounderPerc(jj) ;
        
        % add to the table the percentile
        bmiData.(confVar)(locPercentile) = percVar ;
        
        % increase by 10
        percVar = percVar + 5 ;
        
    end
    
    % convert to a categorical array
    bmiData.(confVar) = categorical(bmiData.(confVar)) ;
    
    % make two separate data set
    bmiDataBlack = bmiData(bmiData.Race == 'AFR', :) ;
    bmiDataEUR = bmiData(bmiData.Race == 'EUR', :) ;
    
    % get a random samples of the EUR equal to the black
    s = RandStream('dsfmt19937') ;
    smallSample = randsample(s,height(bmiDataEUR),height(bmiDataBlack));
    bmiDataEUR = bmiDataEUR(smallSample,:);

    % get the group stats
    plotDataBlack = grpstats(bmiDataBlack,confVar,{'mean','sem'}, ...
        'DataVars',bmiShortTitles);
    plotDataEUR = grpstats(bmiDataEUR,confVar,{'mean','sem'}, ...
        'DataVars',bmiShortTitles);
    
    % plot the variables in a loop
    for ii = 1:length(bmiShortTitles)
        
        % plot the figure for both AFR and EUR
        figure()
        errorbar(plotDataBlack.(confVar), ...
            plotDataBlack.(['mean_',bmiShortTitles{ii}]), ...
            plotDataBlack.(['sem_', bmiShortTitles{ii}]),'-o',...
            'MarkerSize',10,...
            'MarkerEdgeColor',groupColors(1,:),...
            'MarkerFaceColor',groupColors(1,:), ...
            'LineWidth',2)
        
        hold on
        
        errorbar(plotDataEUR.(confVar), ...
            plotDataEUR.(['mean_',bmiShortTitles{ii}]), ...
            plotDataEUR.(['sem_', bmiShortTitles{ii}]),'-o',...
            'MarkerSize',10,...
            'MarkerEdgeColor',groupColors(2,:),...
            'MarkerFaceColor',groupColors(2,:),...
            'LineWidth',2)
        
        % annotate the plot
        set(gca,'FontSize',14,'LineWidth',1.5,'Box','off')
        xlabel([curConf,' Percentile'])
        ylabel([bmiShortTitles{ii}, ' (', bmiVarUnits{ii}, ')'])
        legend({'AFR','EUR'},'Location','SouthWest')
        
        % change the postion of the legend from height 
        if strcmp(curConf,'Height')
           legend({'AFR','EUR'},'Location','SouthEast')
        end
        
        % add the title
        title(bmiFigTitles{ii},'FontSize',16,'FontWeight','bold')
        
        hold off
        
        % save the figure
        if toSave == true
            % create Figure Name to filenames
            name = [ replace(bmiFigTitles{ii},'/','-'), ' ', curConf] ;
            
            % set the properties of the figure and save it
            set(gcf,'paperunits','centimeters','paperposition',[0 0 30 10])
            print(name,'-dpng','-r300','-painters') % -r600
        end
        
    end
end

% ***************** also add some scatter plots *******************
figure()
gscatter(bmiData.FVC,bmiData.BMI,bmiData.Race,groupColors,'..',10)
set(gca,'LineWidth',1,'Box','off','FontWeight','bold')
xlabel('FVC') ;
ylabel('BMI') ;
title('Correlation Between FVC and BMI','FontSize',14','FontWeight','bold')
hold on

% add a reference line and the correlation coefficent
addReferenceLineToPlot(bmiData.FVC,bmiData.BMI)
legend({'Blakcs','EUR'} ,'Location','Best')

figure()
gscatter(bmiData.FVC,bmiData.Age,bmiData.Race,groupColors,'..',10)
set(gca,'LineWidth',1,'Box','off','FontWeight','bold')
xlabel('FVC') ;
ylabel('Age') ;
title('Correlation Between FVC and Age','FontSize',14','FontWeight','bold')
hold on

% add a reference line and the correlation coefficent
addReferenceLineToPlot(bmiData.FVC,bmiData.Age)
legend({'Blakcs','EUR'} ,'Location','Best')

% remove some of the variables 
clear ii locPercentile bmiPerc percVar plotDataBlack plotDataEUR ...
    smallSample ageDataEUR ageDataBlack curConf kk jj bmiShortTitles ...
    bmiVarUnits bmiFigTitles s

%% Clinical Data to Table to Plot the Co-variates and Demographics
% 
% 'WheezeOrWhistlingInTheChestInLastYear'
% 
% 'ShortnessOfBreathWalkingOnLevelGround'
% 
% 'ChestPainOrDiscomfortWalkingNormally'
% 
% 'ChestPainFeltOutsidePhysicalActivity'
% 
% 'ChestPainOrDiscomfortWhenWalkingUphillOrHurrying'

smallData = data(:,[{'EthnicGroupBW','WheezeOrWhistlingInTheChestInLastYear',...
    'ShortnessOfBreathWalkingOnLevelGround','ChestPainOrDiscomfortWalkingNormally',...
    'ChestPainFeltOutsidePhysicalActivity',...
    'ChestPainOrDiscomfortWhenWalkingUphillOrHurrying',...
    'BloodClot_DVT_Bronchitis_Emphysema_Asthma_Rhinitis_Eczema_Aller'}, ...
    myVars]) ;

for ii = 1:width(smallData)
    if iscell(smallData.(ii))
        smallData.(ii) = categorical(smallData.(ii));
    end
end

% first create a glm model comparing EUR and black for FVC
modelAge = ['ForcedVitalCapacity_FVC_',...
    '~EthnicGroupBW',...
    '+WheezeOrWhistlingInTheChestInLastYear', ...
    '+ShortnessOfBreathWalkingOnLevelGround',...
    '+ChestPainOrDiscomfortWhenWalkingUphillOrHurrying',...
    '+BloodClot_DVT_Bronchitis_Emphysema_Asthma_Rhinitis_Eczema_Aller'] ;
glmAll = fitglm(data,modelAge,'Distribution',"normal");
glmBoth = glmAll.Coefficients 

% save the results to excel
writetable(glmBoth, 'Logistic Regression Results.xlsx','Sheet','EUR and AFR')


% create a name for the respiratory variables
respVars =  {'WheezeOrWhistlingInTheChestInLastYear', ...
    'ShortnessOfBreathWalkingOnLevelGround',...
'ChestPainOrDiscomfortWalkingNormally','ChestPainFeltOutsidePhysicalActivity',...
'ChestPainOrDiscomfortWhenWalkingUphillOrHurrying'} ;

% run the text between for ethinicGroup, the categorical features of
% interest and the lung measurement data return only the yes and no answers
curData = smallData ;
for ii = 1 % length(respVars)
    curData = curData(ismember(curData.(respVars{ii}),{'No','Yes'}), :) ;
    curData.(respVars{ii}) = categorical(cellstr(curData.(respVars{ii}))) ;
end

% get only the AFR and data of EUR
logiDataAFR = curData(curData.EthnicGroupBW == "AFR",:) ;
logiDataEUR = curData(curData.EthnicGroupBW == "EUR",:) ;

% specify the model specs
modelAge = ['ForcedVitalCapacity_FVC_',...
    '~WheezeOrWhistlingInTheChestInLastYear', ...
    '+ShortnessOfBreathWalkingOnLevelGround',...
    '+ChestPainOrDiscomfortWhenWalkingUphillOrHurrying'] ;

glmAFR = fitglm(logiDataAFR,modelAge,'Distribution',"normal") ;
coeffAFR = glmAFR.Coefficients ;
glmEUR = fitglm(logiDataEUR,modelAge,'Distribution',"normal") ;
coeffEUR = glmEUR.Coefficients ;

% save the data to excel
writetable(coeffEUR, 'Logistic Regression Results.xlsx','Sheet','EUR')
writetable(coeffAFR, 'Logistic Regression Results.xlsx','Sheet','AFR')
writetable(smallData,'Respiratory Conditions Tableau.xlsx')

clear logiDataAFR logiDataEUR logDataEUR logDataAFR jj modelAge ...
    results timeElapsed uk X Y Z coeffEUR coeffAFR curData ans

%% ********************** Reviewer Comments **************************

% Please discuss how the 5 year difference in age might have affected these
% findings. 5 years could be significant in terms of lung health – younger
% AFR individuals might mean less variability in the lung function
% measurements, which could have affected power to detect associations as
% well.

% first create a glm model comparing EUR and black for FVC
modelAge = ['ForcedVitalCapacity_FVC_',...
    '~EthnicGroupBW','+AgeAtRecruitment'] ;
glmAll = fitglm(data,modelAge,'Distribution',"normal");
glmBoth = glmAll.Coefficients 

figure()
plotSlice(glmAll)

figure()
plotDiagnostics(glmAll)

figure()
plotResiduals(glmAll)

%% ***************** Here is where I run the model **********************

% Specify the model using a formula that allows up to two-way interactions
% between the variables age, weight, and sex. Smoker is the response
% variable. ~EthnicGroupBW','+AgeAtRecruitment
myVarsLong = {'ForcedExpiratoryVolumeIn1_second_FEV1_',...
    'ForcedVitalCapacity_FVC_',...
    'PeakExpiratoryFlow_PEF_'} ;

% get the glm data 
glmData = data(:,[{'EthnicGroupBW','AgeAtRecruitment','StandingHeight'},...
    myVarsLong]) ;
glmData.Properties.VariableNames(1) = "Ethnicity";

% let us change the variables 
for ii = 4:6
    glmData.Properties.VariableNames(ii) = myVars(ii-3) ;
end

% run the generalised model hares
for ii = 1:length(myVars)

    % here are the model specs 
    modelAge = sprintf('%s ~ Ethnicity*AgeAtRecruitment',myVars{ii}) ;
    modelHeight = sprintf('%s ~ Ethnicity*StandingHeight',myVars{ii});

    % Fit a generalised linear model for Age at recruitment
    glmModel.(shortTitles{ii}).age =  ...
        fitglm(glmData,modelAge,'Distribution','normal') ; 

    % Fit a generalised linear model for Age at recruitment
    glmModel.(shortTitles{ii}).height = ...
        fitglm(glmData,modelHeight,'Distribution','normal') ;
end

clear myVarsLong modelAge modelHeight ii

%% Correlation Between the variable of interest

figure()
gscatter( smallData.ForcedVitalCapacity_FVC_ , ...
    smallData.ForcedExpiratoryVolumeIn1_second_FEV1_, ...
    smallData.EthnicGroupBW , groupColors, '..',10)
set(gca,'LineWidth',1,'Box','off','FontWeight','bold')
xlabel('FVC') ;
ylabel('FEV1') ;
title('Correlation Between FVC and FEV1','FontSize',14','FontWeight','bold')
hold on

% add a reference line and the correlation coefficent
addReferenceLineToPlot(smallData.ForcedVitalCapacity_FVC_ , ...
    smallData.ForcedExpiratoryVolumeIn1_second_FEV1_ )
legend({'Blakcs','EUR'} ,'Location','Best')

figure()
gscatter( smallData.ForcedVitalCapacity_FVC_,...
    smallData.PeakExpiratoryFlow_PEF_,...
    smallData.EthnicGroupBW , groupColors, '..',10)
set(gca,'LineWidth',1,'Box','off','FontWeight','bold')
xlabel('FVC') ;
ylabel('PEF') ;
title('Correlation Between FVC and PEF','FontSize',14','FontWeight','bold')
hold on

% add a reference line and the correlation coefficent
addReferenceLineToPlot(smallData.ForcedVitalCapacity_FVC_ , ...
    smallData.PeakExpiratoryFlow_PEF_)
legend({'Blakcs','EUR'} ,'Location','Best')

%% Some Internal Function

function pSuperScript = convertPValue2SuperScript(p)
    % converts the number to scientific superscript for printing on a
    % figure
    pS = num2str(p) ;

    % get the first number
    firstNumbers = extractBefore(pS,'e') ;

    % check if there is a decimal place. then only get the first 4 numbers
    if contains( firstNumbers  ,'.')
        firstNumbers = firstNumbers(1:4) ;
    end

    % get the correctly formated p value
    pSuperScript = sprintf('= %s x 10^{%d}', firstNumbers, ...
        str2double(extractAfter(pS, 'e') )) ;

    % if the p value is large
    if p > 0.0001
        pSuperScript = sprintf('= %0.4f', p) ;
    elseif p == 0
        pSuperScript = sprintf('< 1 x 10^{%d}', -350) ;
    end

end


function addReferenceLineToPlot(xVariable, yVariable)
% This function add a reference line to the plot and the R square value

% replace the missing values with mean
xVariable = fillmissing(xVariable,"linear");
yVariable = fillmissing(yVariable,"linear");

% preform a regression calculation and add it to the plot
X = [ones(length(xVariable),1) xVariable] ;
b1 = X\yVariable    ;     yCalc1 = X*b1 ;
plot(xVariable,yCalc1,'LineWidth',1.5,'Color','k')

% calculate the pearson's linear correation  and add it to the plot
[r2 , pR ] = corr(xVariable,yVariable, 'Type','Pearson',...
    'Rows','complete');
if pR < 0.0001
    text(0.4, 0.9 , strcat( sprintf( ...
        "R = %0.2f, P ", r2), convertPValue2SuperScript(pR)  ), ...
        'Units','normalized','FontSize',11,'FontWeight','bold')
else
    text(0.4, 0.9 ,sprintf('R = %0.2f, P = %0.4f', r2, pR), ...
        'Units','normalized','FontSize',11,'FontWeight','bold')
end
end

%% ====================== Internal Functions =================

function createLegendInternal(yPoint, xStart, legendLabels , plotColors,...
    myLgdTitle , fontSizes ,rectAndTextBox)

% specificy the y values starts and mode of actions for the drugs
% yPoint = 0.830 ; xStart = 0.1023 ;
xStartText = xStart + 0.01 ;
yPointTitle = yPoint + 0.03 ;

% specify the font size to be used in the plot
if ~exist('fontSizes','var')
    fontSizes = [10, 12] ;
end

% specifiy the rectangle and text length
if ~exist('rectAndTextBox','var')
    rectAndTextBox = [0.018 ,0.12] ;
end

% check for errors
if ~isnumeric(yPoint) || ~isnumeric(xStart)
    error('Both yPoint and xStarts should be numeric values')
elseif yPoint > 1 || xStart > 1
    error('Both yPoint and xStarts should be less than 1')
elseif ~isnumeric(plotColors)
    error('plot Color should be numeric')
end

if size(plotColors,1) ~= length(legendLabels)
    error('There should be a color for each legend names')
end

if iscategorical( legendLabels)
    legendLabels = categories(legendLabels);
end

for ii = 1:length(legendLabels)
    % add the legend color
    annotation('rectangle',[xStart yPoint rectAndTextBox(1) 0.023],...
        'EdgeColor', plotColors(ii,:), ...
        'FaceColor', plotColors(ii,:));
    
    % add the legend text
    annotation('textbox',[xStartText yPoint rectAndTextBox(2) 0.0230],...
        'String',legendLabels{ii},'FontSize',fontSizes(1),...
        'FontName','Helvetica Neue','FitBoxToText','off',...
        'EdgeColor',[1 1 1],'BackgroundColor',[1 1 1] , ...
        'VerticalAlignment','middle','FontWeight','normal')
    
    % move the y point down
    yPoint = yPoint - 0.03 ;
end

% add the title
annotation('textbox',[xStart yPointTitle rectAndTextBox(2) 0.0230],...
    'String', myLgdTitle,'FontSize',fontSizes(2),...
    'FontName','Helvetica Neue','FitBoxToText','off',...
    'EdgeColor',[1 1 1],'BackgroundColor',[1 1 1] , ...
    'VerticalAlignment','middle','FontWeight','bold',...
    'HorizontalAlignment','left');

end

% **************
%% 

%% The Internal Function 

function ManhattanPlotGrey( filename, varargin )

%   This function takes a GWAS output file  and plots a Manhattan Plot.

% ARGUMENTS:
%   sex: defaults to 0, set to 1 to include sex chromosomes

%   sig: defaults to 5e-8, significance threshold for horizontal line

%   vert: defaults to 0, set to 1 for tick labels have chr in them and are
%   written upwards

%   labels: defaults to [-1,-1]. Set to [x,y] to label the top SNP on each
%   locus with p<x. Locus defined by windows of y base pairs.

%   outfile: defaults to filename of input file, set to something else to
%   change name of output file

%   title: defaults to 'Manhattan Plot'. Fairly self-explanatory

%   save: defaults to 1, set to 0 to disable saving to save time

%   format: defaults to PLINK, options {PLINK,BOLT-LMM,SAIGE} This is used
%   to identify the correct column headings. If using anything else, rename
%   the header line so that CHR, BP, and P are the headers for chromosome,
%   base pair and p value, and the default PLINK should catch it

%   Usage:
%   ManhattanPlot('gwas.assoc.fisher',varargin) will take the association
%   analysis in 'gwas.assoc.fisher, generate a Manhattan Plot, and store it
%   in gwas_ManhattanPlot.png, which should be publication-ready, and
%   gwas_ManhattanPlot.fig for minor readjustments in MATLAB's GUI. This is
%   fine as .fisher is a PLINK format file. For a SAIGE output (and
%   significance threshold 5e-9) with sex chromosomes shown, and the top
%   hit for each SNP within 1 mb and below p=1e-6 labelled, use
%   ManhattanPlot('gwas.saige','format','SAIGE','sig',5e-9,'sex',1,
%   'labels',[1e-6,1000000])

%   Tested using assoc, assoc.logistic and assoc.fisher files generated by
%   Plink 1.7. Also tested using BOLT-LMM and SAIGE.
%
%   Harry Green (2019) Genetics of Complex Traits, University of Exeter

% Reading in optional arguments and input file

p = inputParser;
defaultsex = 0;
defaultvert = 0;
defaultsig = 5e-8;
defaultsave = 1;
defaultoutfile = filename;
defaulttitle= 'Manhattan Plot';
defaultlabels = [-1, -1];
defaultformat = 'PLINK';
expectedformat = {'PLINK','BOLT-LMM','SAIGE'};

% here is the number of expected plots 
defaultplotNumber = 1;
defaultplotCounter = 1 ;

addRequired(p,'filename');
% addOptional(p,'outfile',defaultoutfile,@ischar);
addOptional(p,'sex',defaultsex,@isnumeric);
addOptional(p,'vert',defaultvert,@isnumeric);
addOptional(p,'labels',defaultlabels,@isnumeric);
addOptional(p,'sig',defaultsig,@isnumeric);
addOptional(p,'save',defaultsave,@isnumeric);
addParameter(p,'format',defaultformat,...
    @(x) any(validatestring(x,expectedformat)));
addParameter(p,'outfile',defaultoutfile,@ischar);
addParameter(p,'title',defaulttitle,@ischar);
addParameter(p,'plotColor',@isnumeric);

addParameter(p,'plotNumber',defaultplotNumber,@isnumeric);
addParameter(p,'plotCounter',defaultplotCounter,@isnumeric);

parse(p,filename,varargin{:});

filename=p.Results.filename;
sex=p.Results.sex;
vert=p.Results.vert;
sig=p.Results.sig;
save=p.Results.save;
format=p.Results.format;
outfile=p.Results.outfile;
labels=p.Results.labels;
plottitle=p.Results.title;
plotColor=p.Results.plotColor ;

plotNumber=p.Results.plotNumber ;
plotCounter=p.Results.plotCounter ;

% ################ Check if the input is a table for a file name #########

if ischar(filename)
    % MATLAB refuses to read file formats like .assoc, so first rename as
    % .txt
    copyfile(filename, strcat(filename,'.txt'));
    opts = detectImportOptions(strcat(filename,'.txt'),'NumHeaderLines',0);
    
    % the columns of T will be the headers of the assoc file
    T = readtable(strcat(filename,'.txt'),opts);
    
    %delete unwanted .txt file that we didn't want anyway
    delete(strcat(filename,'.txt'))
elseif istable(filename)
    T = filename ;
end


% File format defintitions
% This section uses the default column headings from PLINK, BOLT-LMM and
% SAIGE. Easy to modify for other software packages

% check the CHR variable exist then change to PLINK format
if ~ismember(T.Properties.VariableNames,'CHR')
    T.Properties.VariableNames([1,2]) = {'CHR','BP'} ;
end

if strcmp(format,'PLINK')
    % These are the only useful columns for a manhattan plot, chromosome,
    % base pair, p value
    tab2=[T.CHR,T.BP,T.P]; 
end

if strcmp(format,'BOLT-LMM')
    % These are the only useful columns for a manhattan plot, chromosome,
    % base pair, p value
    tab2=[T.CHR,T.BP,T.P_BOLT_LMM]; 
end

if strcmp(format,'SAIGE')
    % These are the only useful columns for a manhattan plot, chromosome,
    % base pair, p value
    tab2=[T.CHR,T.POS,T.p_value]; 
end

%

if sex==0
    tab2=tab2(tab2(:,1)<23,:); %remove chr23 if not wanting sex chromosomes
end

% variable to track base pairs, this helps the plotting function to know
% where to start plotting the next chromosome
bptrack=0; 
tab2(tab2(:,3)==0,:)=[];

% get the signficant loci based on the pvalues
tabSig = tab2(tab2(:,3) < sig, :) ;

% get the  non significant loci
tab2 = tab2(tab2(:,3) >= sig, :) ;

% % get the plot
% figure(plotNumber)

% produce the plot in a loop
for i=1:22+sex
    hold on
    
    %a scatterplot. On the x axis, the base pair number + bptrack, which
    %starts at 0. On the y axis, - log 10 p
    plot( tab2(tab2(:,1)==i,2)+max(bptrack),...
        -log10(tab2(tab2(:,1)==i,3)), '.' , ...
        'MarkerEdgeColor',[0.5 0.5 0.5],...
        'MarkerFaceColor',[0.5 0.5 0.5],...
        'MarkerSize',10);
    
    % produce a scatter plot of the significant variants
    plot( tabSig(tabSig(:,1)==i,2)+max(bptrack),...
        -log10(tabSig(tabSig(:,1)==i,3)), '.' , ...
        'MarkerEdgeColor',plotColor,...
        'MarkerFaceColor',plotColor,...
        'MarkerSize',10);
    
    %this updates bptrack by adding the highest base pair number of the
    %chromosome. At the end, we should get the total number of base pairs
    %in the human genome. All values of bptrack are stored in a vector.
    %They're useful later
    bptrack=[bptrack,max(bptrack)+max(max(tab2(tab2(:,1)==i,:)))]; 
end

if plotNumber > plotCounter
    return
end

%if strongest hit is 1e-60, plot window goes up to 1e-61
if isempty(tabSig)
    %and down to the highest p value.
    ylimit=max(max(-log10(tab2(:,3)))+1);
    ylimitmin=min(min(-log10(tab2(:,3))));
else
    %and down to the highest p value.
    ylimit=max(max(-log10(tabSig(:,3)))+1);
    ylimitmin=min(min(-log10(tabSig(:,3))));
end

%genome wide significant line, uses the sig optional argument which
%defaults to 5e-8
plot([0,max(bptrack)],-log10([sig,sig]),'k--') 

% ylabel('$-\log_{10}(p)$','Interpreter','latex')
ylabel('-log_1_0(p)')

% xlim([0,max(bptrack)]) % I HAVE CHANGED THIS 

% find the previous ylimit
curYlim = get(gca, 'YLim');

% compared with the new yLim
if max(ylimit,-log10(sig/5)) > curYlim
    %y axis will always go to the significance threshold/5 at least
    ylim([0,max(ylimit,-log10(sig/5))]) % floor(ylimitmin)
end

%this calculates the moving average of bptrack and uses them as chromosome
%label markers. This puts them in the middle of each chromosome.
M=movmean(bptrack,2); 
xticks(M(2:end));

% Rotation section
% this section of the code changes the x axis label orientation.

if vert ==0
    xlabel('Chromosome')% ,'Interpreter','latex')
    xticklabels( 1:23 );
end
if vert ==1
    xtickangle(90)
    xticklabels( {'chr1','chr2','chr3','chr4','chr5','chr6','chr7',...
        'chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr15',...
        'chr16','chr17','chr18','chr19','chr20','chr21','chr22','chrXY'});
end

% Annotation section
% this section of the code annotates the top SNP on each locus

% Neither of these should be negative. Defaults to -1 -1
if ~(labels(1)==-1||labels(2)==-1) 
    labellim=labels(1);
    labelprox=labels(2);
    
    %as nothing else will be used for labelling, the reduced table is
    %easier to manage
    tab3=tab2(tab2(:,3)<labellim,:); 
    
    %easier if they're in ascending order of p value;
    [~,index]=sort(tab3(:,3));
    tab3=tab3(index,:);
    
    %this becomes a list of SNPs that have been labelled. It starts with a
    %dummy SNP to give the first SNP something to compare with. This will
    %always fail so the first label always gets added
    labelledsnps=[0,0,0]; 
    
    for i=1:size(tab3,1)
        test=abs(tab3(i,:)-labelledsnps);
        %if there are no snps on the same chromosome and within 1 MB
        if sum(test(:,1)==0&test(:,2)<labelprox)==0 
            
            %this plots a square over the corresponding dot
            plot(bptrack(tab3(i,1))+tab3(i,2),-log10(tab3(i,3)),'ks',...
                'MarkerSize',8,'LineWidth',1) 
            
            %this puts together a strong of chrx:y
            labels=strcat('chr',string(tab3(i,1)),':',...
                string(bptrack(tab3(i,1))+tab3(i,2))); 
            
            %this plots the label above and to the left of the dot, shifted
            %up by 0.05 to give some space
            text(bptrack(tab3(i,1))+tab3(i,2),-log10(tab3(i,3))+0.05,...
                char(labels),'VerticalAlignment','bottom',...
                'HorizontalAlignment','left') %,'Interpreter','latex') 
            labelledsnps=[labelledsnps;tab3(i,:)];
        end
    end
end
% finalising file

% takes the output file name up to the first decimal point
a = strsplit(outfile,'.'); 

% my code 
set(gca,'FontSize',12,'LineWidth',0.5)

box on
title(plottitle,'fontsize',14)
% title(plottitle,'Interpreter','latex','fontsize',14)

if save~=0
    name=[a{1},'_ManhattanPlot.fig']; %adds _ManhattanPlot to filenames
    savefig(name);
    set(gcf, 'paperunits', 'centimeters', 'paperposition', [0 0 30 10])
    print([a{1},'_ManhattanPlot.png'],'-dpng','-r300', '-vector') % -r600
end
end

%% Another internal function 

function varargout = venn (varargin)

%VENN   Plot 2- or 3- circle area-proportional Venn diagram
%
%  venn(A, I)
%  venn(Z)
%  venn(..., F)
%  venn(..., 'ErrMinMode', MODE)
%  H = venn(...)
%  [H, S] = venn(...)
%  [H, S] = venn(..., 'Plot', 'off')
%  S = venn(..., 'Plot', 'off')
%  [...] = venn(..., P1, V1, P2, V2, ...) 
%
%venn(A, I) by itself plots circles with total areas A, and intersection
%area(s) I. For two-circle venn diagrams, A is a two element vector of circle 
%areas [c1 c2] and I is a scalar specifying the area of intersection between 
%them. For three-circle venn diagrams, A is a three element vector [c1 c2 c3], 
%and I is a four element vector [i12 i13 i23 i123], specifiying the 
%two-circle intersection areas i12, i13, i23, and the three-circle
%intersection i123.
%
%venn(Z) plots a Venn diagram with zone areas specified by the vector Z. 
%For a 2-circle venn diagram, Z is a three element vector [z1 z2 z12]
%For a 3-circle venn, Z is a 7 element vector [z1 z2 z3 z12 z13 z23 z123]
%
%venn(..., F) specifies optional optimization options. VENN uses FMINBND to
%locate optimum pair-wise circle distances, and FMINSEARCH to optimize
%overall three-circle alignment. F is a structure with fields specifying
%optimization options for these functions. F may be a two-element array of
%structures, in which case the first structure is used for FMINBND
%function calls, and the second structure is used for FMINSEARCH function
%calls.
%
%venn(..., 'ErrMinMode', MODE)
%Used for 3-circle venn diagrams only. MODE can be 'TotalError' (default), 
%'None', or 'ChowRodgers'. When ErrMinMode is 'None', the positions and 
%sizes of the three circles are fixed by their pairwise-intersections, 
%which means there may be a large amount of error in the area of the three-
%circle intersection. Specifying ErrMinMode as 'TotalError' attempts to 
%minimize the total error in all four intersection zones. The area of the 
%three circles are kept constant in proportion to their populations. The 
%'ChowRodgers' mode uses the the method proposed by Chow and Rodgers 
%[Ref. 1] to draw 'nice' three-circle venn diagrams which appear more 
%visually representative of the desired areas, although the actual areas of 
%the circles are allowed to deviate from requested values.
%
%H = venn(...) returns a two- or three- element vector to the patches 
%representing the circles. 
%
%[H, S] = venn(...) returns a structure containing descriptive values
%computed for the requested venn diagram. S is a structure with the
%following fields, where C is the number of circles (N = 2 or 3), Z is
%the number of zones (Z = 3 or 7), and I is the number of intersection 
%areas (1 or 4)
%
% Radius            C-element vector of circle radii
%
% Position          C*2 array of circle centers
%
% ZoneCentroid      Z*2 array of zone centroids (Can be used for labeling)
%
% CirclePop         C-element vector of supplied circle populations. 
%                   (I.e., the 'true' circle areas)
%
% CircleArea        C-element of actual circle areas
%
% CircleAreaError   = (CircleArea-CirclePop)/CirclePop
%
% IntersectPop      I-element vector of supplied intersection populations
%                   (I.e., the 'true' intersection areas)
%
% IntersectArea     I-element vector of actual intersection areas
%
% IntersectError    = (IntersectArea-IntersectPop)/IntersectPop
%
% ZonePop           Z-element vector of supplied zone populations. (I.e.
%                   'true' zone areas
%
% ZoneArea          Z-element vector of actual zone areas.
%
% ZoneAreaError     = (ZoneArea-ZonePop)/ZonePop
% 
%
%[H, S] = venn(..., 'Plot', 'off')
%S = venn(..., 'Plot', 'off')

% Returns a structure of computed values, without plotting the diagram.
% This which can be useful when S is used to draw custom venn diagrams or
% for exporting venn diagram data to another application. When Plot is set
% to off, the handles vector H is returned as an empty array.
% Alternatively, the command

% S = venn(..., 'Plot', 'off) will return only the output structure.

%[...] = venn(..., P1, V1, P2, V2, ...) 

%Specifies additional patch settings in standard Matlab parameter/value
%pair syntax. Parameters can be any valid patch parameter. Values for patch
%parameters can either be single values, or a cell array of length
%LENGTH(A), in which case each value in the cell array is applied to the
%corresponding circle in A.

%Examples
%
%   %Plot a simple 2-circle venn diagram with custom patch properties
%   figure, axis equal, axis off
%   A = [300 200]; I = 150;
%   venn(A,I,'FaceColor',{'r','y'},'FaceAlpha',{1,0.6},'EdgeColor','black')
%
%   %Compare ErrMinModes
%   A = [350 300 275]; I = [100 80 60 40];
%   figure
%   subplot(1,3,1), h1 = venn(A,I,'ErrMinMode','None');
%   axis image,  title ('No 3-Circle Error Minimization')
%   subplot(1,3,2), h2 = venn(A,I,'ErrMinMode','TotalError');
%   axis image,  title ('Total Error Mode')
%   subplot(1,3,3), h3 = venn(A,I,'ErrMinMode','ChowRodgers');
%   axis image, title ('Chow-Rodgers Mode')
%   set([h1 h2], 'FaceAlpha', 0.6)
%
%   %Using the same areas as above, display the error optimization at each 
%   iteration. Get the output structure.
%   F = struct('Display', 'iter');
%   [H,S] = venn(A,I,F,'ErrMinMode','ChowRodgers','FaceAlpha', 0.6);
%
%   %Now label each zone 
%   for i = 1:7
%       text(S.ZoneCentroid(i,1), S.ZoneCentroid(i,2), ['Zone ' num2str(i)])
%   end
%
%See also patch, bar, optimset, fminbdn, fminsearch
%
%Copyright (C) 2008 Darik Gamble, University of Waterloo.
%dgamble@engmail.uwaterloo.ca
%
%References
%1. S Chow and P Rodgers. Extended Abstract: Constructing Area-Proportional
%   Venn and Euler Diagrams with Three Circles. Presented at Euler Diagrams 
%   Workshop 2005. Paris. Available online: 
%   http://www.cs.kent.ac.uk/pubs/2005/2354/content.pdf
%
%2. S Chow and F Ruskey. Drawing Area-Proportional Venn and Euler Diagrams. 
%   Lecture Notes in Computer Science. 2004. 2912: 466-477. Springer-Verlag. 
%   Available online: http://www.springerlink.com/content/rxhtlmqav45gc84q/
%
%3. MP Fewell. Area of Common Overlap of Three Circles. Australian Government 
%   Department of Defence. Defence Technology and Science Organisation. 2006. 
%   DSTO-TN-0722. Available online:
%   http://dspace.dsto.defence.gov.au/dspace/bitstream/1947/4551/4/DSTO-TN-0722.PR.pdf


%Variable overview
%   A0, A   Desired and actual circle areas
%               A = [A1 A2] or [A1 A2 A3]
%   I0, I   Desired and actual intersection areas
%               I = I12 or [I12 I13 I23 I123]
%   Z0, Z   Desired and actual zone areas
%               Z = [Z1 Z2 Z12] or [Z1 Z2 Z3 Z12 Z13 Z23 Z123]
%   x, y    Circle centers
%               x = [x1 x2] or [x1 x2 x3]
%   r       Circle radii
%               r = [r1 r2] or [r1 r2 r3]
%   d       Pair-wise distances between circles
%               d = d12 or [d12 d13 d23]
    


    %Parse input arguments and preallocate settings
    [A0, I0, Z0, nCirc, fminOpts, vennOpts, patchOpts] = ...
        parseArgsIn(varargin);
    [d, x, y, A, I, Z] = preallocVectors (nCirc);
    zoneCentroids = []; %Will only be calculated if needed
    
    %Circle Radii
    r = sqrt(A0/pi);

    %Determine distance between first circle pair    
    d(1) = circPairDist(r(1), r(2), I0(1), fminOpts(1));
    
    %Position of second circle is now known
    x(2) = d(1); 
    
    %First intersection area 
    I(1) = areaIntersect2Circ(r(1), r(2), d(1));
    
    if nCirc==3
        %Pairwise distances for remaining pairs 1&3 and 2&3
        d(2) = circPairDist(r(1), r(3), I0(2), fminOpts(1)); %d13
        d(3) = circPairDist(r(2), r(3), I0(3), fminOpts(1)); %d23

        %Check triangle inequality
        srtD = sort(d);
        if ~(srtD(end)<(srtD(1)+srtD(2)))
            error('venn:triangleInequality', ...
                'Triangle inequality not satisfied')
        end

        %Guess the initial position of the third circle using the law of
        %cosines
        alpha = acos( (d(1)^2 + d(2)^2 - d(3)^2)  / (2 * d(1) * d(2)) );
        x(3) = d(2)*cos(alpha);
        y(3) = d(2)*sin(alpha);

        %Each pair-wise intersection fixes the distance between each pair
        %of circles, so technically there are no degrees of freedom left in
        %which to adjust the three-circle intersection. We can either try
        %moving the third circle around to minimize the total error, or
        %apply Chow-Rodgers 
        
        switch vennOpts.ErrMinMode
            case 'TotalError'
                %Minimize total intersection area error by moving the third
                %circle
                pos = fminsearch(@threeCircleAreaError, [x(3) y(3)],...
                    fminOpts(2));
                x(3) = pos(1);
                y(3) = pos(2);
            case 'ChowRodgers'
                %note that doChowRodgersSearch updates x and y in this
                %workspace as a nested fcn
                doChowRodgersSearch;
        end

        %Make sure everything is 'up to date' after optimization
        update3CircleData;
        
    end
    
    % Are we supposed to plot?
    if vennOpts.Plot
        if isempty(vennOpts.Parent)
            vennOpts.Parent = gca;
        end
        hVenn = drawCircles(vennOpts.Parent, x, y, r, ...
            patchOpts.Parameters, patchOpts.Values);
    else
        hVenn = [];
    end
    
    %Only determine zone centroids if they're needed 
    %Needed for output structure 
    nOut = nargout;
    if (nOut==1 && ~vennOpts.Plot) || nOut==2
        if nCirc == 2
            %Need to calculate new areas
            A = A0; %Areas never change for 2-circle venn
            Z = calcZoneAreas(2, A, I);
            zoneCentroids = zoneCentroids2(d, r, Z);
        else
            zoneCentroids = zoneCentroids3(x, y, d, r, Z);
        end
    end
        
    %Figure out output arguments
    if nOut==1
        if vennOpts.Plot
            varargout{1} = hVenn;
        else
            varargout{1} = getOutputStruct;
        end
    elseif nOut==2
        varargout{1} = hVenn;
        varargout{2} = getOutputStruct;
    end
        
    
    
    function err = threeCircleAreaError (pos)
        
        x3 = pos(1);
        y3 = pos(2);
        
        %Calculate distances
        d(2) = sqrt(x3^2 + y3^2); %d13
        d(3) = sqrt((x3-d(1))^2 + y3^2); %d23
        
        %Calculate intersections
        %Note: we're only moving the third circle, so I12 is not changing
        I(2:3) = areaIntersect2Circ (r(1:2), r([3 3]), d(2:3)); %I13 and I23
        I(4) = areaIntersect3Circ (r, d); %I123
        
        %Replace 0 (no intersection) with infinite error
        I(I==0) = Inf;
        
        %Error
        err = sum(abs((I-I0)./I0));
        
    end


    function doChowRodgersSearch
        
        %Adapted from Ref. [1]
        
        %Initialize an index matrix to select all 7choose2 zone pairs (21 pairs)
        idx = nchoosek(1:7, 2);
                
        %Which zone-zone pairs are considered equal?
        %Zones within 10% of each other considered equal
        zonePairAreas0 = Z0(idx);
        
        %Percent difference in population between the two members of a pair
        ar0 = 2*abs(zonePairAreas0(:,1)-zonePairAreas0(:,2))./sum(zonePairAreas0, 2)*100;
        eqPairCutoff = 10;  
        pairIsEq = ar0<=eqPairCutoff;
        
        %Calculate allowable range for pairs of zones considered unequal
        if any(~pairIsEq)
            %Sort zone areas
            [zUneqAreas0, zUneqAreasSrtIdx] = sort(zonePairAreas0(~pairIsEq,:), 2);
            
            %Make a real index array out of the inconvenient index sort returns
            n = sum(~pairIsEq);
            zUneqAreasSrtIdx = sub2ind([n,2], [1:n; 1:n]', zUneqAreasSrtIdx);
            
            %rp = (largepopulation/smallpopulation)-1
            rp = zUneqAreas0(:,2)./zUneqAreas0(:,1)-1;
            rpMin = 1 + 0.3*rp;
            rpMax = 1 + 2*rp;
        end
        
        %Preallocate zone error vector
        zoneErr = zeros(1,21); 

        %Initialize independent parameters to search over
        guessParams = [r(1) x(2) r(2) x(3) y(3) r(3)];
        
        %Search!
        pp = fminsearch(@chowRodgersErr, guessParams, fminOpts(2));  
        
        [r(1) x(2) r(2) x(3) y(3) r(3)] = deal(pp(1), pp(2), pp(3), pp(4), pp(5), pp(6));
        
        
        function err = chowRodgersErr (p)
            
            %params = [x2 r2 x3 y3 r3]
            [r(1), x(2), r(2), x(3), y(3), r(3)] = deal(p(1), p(2), p(3), p(4), p(5), p(6));
             
            %After changing x2, r2, x3, y3, and r3, update circle areas,
            %distances, intersection areas, zone areas
            update3CircleData;

            if any(pairIsEq)
                %For zone pairs considered equal, error is equal to square of the
                %distance beyond the cutoff; 0 within cutoff
                zAreas = Z(idx(pairIsEq,:));
                ar = 2*abs(zAreas(:,1)-zAreas(:,2))./sum(zAreas, 2)*100;
                isWithinRange = ar<eqPairCutoff;
                ar(isWithinRange) = 0;
                ar(~isWithinRange) = ar(~isWithinRange) - eqPairCutoff;

                %Amplify error for equal zones with unequal areas
                eqZoneUneqAreaErrorGain = 10;
                ar(~isWithinRange) = ar(~isWithinRange)*eqZoneUneqAreaErrorGain;

                zoneErr(pairIsEq) = ar.^2;
            end

            if any(~pairIsEq)
                %For zone pairs considered unequal, error is equal to square of
                %the distance from the allowable range of rp

                %rp = (largepopulation/smallpopulation)-1
                zUneqPairAreas = Z(idx(~pairIsEq,:));
                
                %Sort based on the population sizes (determined by parent
                %function doChowRodgersSearch)
                zUneqPairAreas = zUneqPairAreas(zUneqAreasSrtIdx);
                rp = zUneqPairAreas(:,2)./zUneqPairAreas(:,1)-1;

                lessThanMin = rp<rpMin;
                moreThanMax = rp>rpMax;
                rp(~lessThanMin & ~moreThanMax) = 0;
                
                %Determine how far out of range errors are
                rp(lessThanMin) = rp(lessThanMin) - rpMin(lessThanMin);
                rp(moreThanMax) = rp(moreThanMax) - rpMax(moreThanMax);          

                %Consider the case where rp < rpMin to be more
                %erroneous than the case where rp > rpMax 
                tooSmallErrorGain = 10;
                rp(lessThanMin) = rp(lessThanMin)*tooSmallErrorGain;

                zoneErr(~pairIsEq) = rp.^2;
            end
            
            %Total error
            err = sum(zoneErr);
            
        end %chowRodgersErr
        
    end %doChowRodgersSearch

    function update3CircleData
        
        %Circle areas
        A = pi*r.^2;

        %Calculate distances
        d(1) = abs(x(2)); %d12
        d(2) = sqrt(x(3)^2 + y(3)^2); %d13
        d(3) = sqrt((x(3)-d(1))^2 + y(3)^2); %d23

        %Calculate actual intersection areas
        I(1:3) = areaIntersect2Circ (r([1 1 2]), r([2 3 3]), d); %I12, I13, I23
        I(4) = areaIntersect3Circ (r, d); %I123

        %Calculate actual zone areas
        Z = calcZoneAreas(3, A, I);

    end

    function S = getOutputStruct
               
        S = struct(...
            'Radius'                ,r                      ,...        
            'Position'              ,[x' y']                ,...
            'ZoneCentroid'          ,zoneCentroids          ,...
            'CirclePop'             ,A0                     ,...
            'CircleArea'            ,A                      ,...
            'CircleAreaError'       ,(A-A0)./A0             ,...
            'IntersectPop'          ,I0                     ,...
            'IntersectArea'         ,I                      ,...
            'IntersectError'        ,(I-I0)./I0             ,...
            'ZonePop'               ,Z0                     ,...
            'ZoneArea'              ,Z                      ,...
            'ZoneAreaError'         ,(Z-Z0)./Z0             );  
        end

end %venn

        
function D = circPairDist (rA, rB, I, opts)
    %Returns an estimate of the distance between two circles with radii rA and
    %rB with area of intersection I
    %opts is a structure of FMINBND search options
    D = fminbnd(@areadiff, 0, rA+rB, opts);
    function dA = areadiff (d)
        intersectArea = areaIntersect2Circ (rA, rB, d);
        dA = abs(I-intersectArea)/I;
    end
end

function hCirc = drawCircles(hParent, xc, yc, r, P, V,c)

    hAx = ancestor(hParent, 'axes');
    nextplot = get(hAx, 'NextPlot');
    
    %P and V are cell arrays of patch parameter/values
    xc = xc(:); yc = yc(:);     %Circle centers
    r = r(:);                   %Radii
    n = length(r);              
    
    %Independent parameter
    dt = 0.05;
    t = 0:dt:2*pi;

    % Origin centered circle coordinates
    X = r*cos(t);
    Y = r*sin(t);
    
    hCirc = zeros(1,n);
    % c = [0.00,0.45,0.74; 0.85,0.33,0.10; 0.47,0.67,0.19];
    c = [0,0.450,0.740;0.640,0.080,0.18] ; %blue and orange-red
    % c = [ 0.85,0.33,0.10; 0.47,0.67,0.19  ]; % green and orange-brown
    
    % c = {'r', 'g', 'b'};                        %default colors
    fa = {0.6, 0.6, 0.6};                       %default face alpha
    tag = {'Circle1', 'Circle2', 'Circle3'}; 	%default tag
    
    for i = 1:n
        xx = X(i,:)+xc(i);  
        yy = Y(i,:)+yc(i);
        
        if iscell(c)
            % if the color is cell array
            hCirc(i) = patch (xx, yy, c{i}, 'FaceAlpha', fa{i},...
                'Parent', hParent, 'Tag', tag{i});
        else % or my colors
            hCirc(i) = patch (xx, yy, c(i,:), 'FaceAlpha', fa{i},...
                'Parent', hParent, 'Tag', tag{i} ,'LineWidth',1);
        end
        
        if i==1
            set(hAx, 'NextPlot', 'add');
        end
    end
    set(hAx, 'NextPlot', nextplot);
    
    % make the axis invisible 
    set(findobj(gcf, 'type','axes'), 'Visible','off')

%     % Custom patch parameter values
%     if ~isempty(P)
% 
%         c = cellfun(@iscell, V);
% 
%         %Scalar parameter values -- apply to all circles
%         if any(~c)
%             set(hCirc, {P{~c}}, {V{~c}});
%         end
% 
%         %Parameters values with one value per circle
%         if any(c)
%             %Make sure all vals are column cell arrays
%             V = cellfun(@(val) (val(:)), V(c), 'UniformOutput', false);
%             set(hCirc, {P{c}}, [V{:}])
%         end
%     end
    
end %plotCircles

 
function A = areaIntersect2Circ (r1, r2, d)
    %Area of Intersection of 2 Circles
    %Taken from [2]
    
    alpha = 2*acos( (d.^2 + r1.^2 - r2.^2)./(2*r1.*d) );
    beta  = 2*acos( (d.^2 + r2.^2 - r1.^2)./(2*r2.*d) );
    
    A =    0.5 * r1.^2 .* (alpha - sin(alpha)) ...  
         + 0.5 * r2.^2 .* (beta - sin(beta));
    
end

function [A, x, y, c, trngArea] = areaIntersect3Circ (r, d)
    %Area of common intersection of three circles
    %This algorithm is taken from [3]. 
    %   Symbol    Meaning
    %     T         theta
    %     p         prime
    %     pp        double prime
        
    %[r1 r2 r3] = deal(r(1), r(2), r(3));
    %[d12 d13 d23] = deal(d(1), d(2), d(3));

    %Intersection points
    [x,y,sinTp,cosTp] = intersect3C (r,d);
    
    if any(isnan(x)) || any(isnan(y))
        A = 0;
        %No three circle intersection
        return
    end
    
    %Step 6. Use the coordinates of the intersection points to calculate the chord lengths c1,
    %c2, c3:
    i1 = [1 1 2];
    i2 = [2 3 3];
    c = sqrt((x(i1)-x(i2)).^2 + (y(i1)-y(i2)).^2)';

    %Step 7: Check whether more than half of circle 3 is included in the circular triangle, so
    %as to choose the correct expression for the area
    lhs = d(2) * sinTp;
    rhs = y(2) + (y(3) - y(2))/(x(3) - x(2))*(d(2)*cosTp - x(2));
    if lhs < rhs
        sign = [-1 -1 1];
    else
        sign = [-1 -1 -1];
    end
    
    %Calculate the area of the three circular segments.
    ca = r.^2.*asin(c/2./r) + sign.*c/4.*sqrt(4*r.^2 - c.^2);

    trngArea = 1/4 * sqrt( (c(1)+c(2)+c(3))*(c(2)+c(3)-c(1))*(c(1)+c(3)-c(2))*(c(1)+c(2)-c(3)) );
    A = trngArea + sum(ca);
    
end

function [x, y, sinTp, cosTp] = intersect3C (r, d)
    %Calculate the points of intersection of three circles
    %Adapted from Ref. [3]
    
    %d = [d12 d13 d23]
    %x = [x12; x13; x23]
    %y = [y12; y13; y23]

    %   Symbol    Meaning
    %     T         theta
    %     p         prime
    %     pp        double prime
    
    x = zeros(3,1);
    y = zeros(3,1);
     
    %Step 1. Check whether circles 1 and 2 intersect by testing d(1)
    if ~( ((r(1)-r(2))<d(1)) && (d(1)<(r(1)+r(2))) )
        %x = NaN; y = NaN;
        %bigfix: no returned values for sinTp, cosTp
        [x, y, sinTp, cosTp] = deal(NaN);
        return
    end

    %Step 2. Calculate the coordinates of the relevant intersection point of circles 1 and 2:
    x(1) = (r(1)^2 - r(2)^2 + d(1)^2)/(2*d(1));
    y(1) = 0.5/d(1) * sqrt( 2*d(1)^2*(r(1)^2 + r(2)^2) - (r(1)^2 - r(2)^2)^2 - d(1)^4 );

    %Step 3. Calculate the values of the sines and cosines of the angles tp and tpp:
    cosTp  =  (d(1)^2 + d(2)^2 - d(3)^2) / (2 * d(1) * d(2));
    cosTpp = -(d(1)^2 + d(3)^2 - d(2)^2) / (2 * d(1) * d(3));
    sinTp  =  (sqrt(1 - cosTp^2));
    sinTpp =  (sqrt(1 - cosTpp^2));

    %Step 4. Check that circle 3 is placed so as to form a circular triangle.
    cond1 = (x(1) - d(2)*cosTp)^2 + (y(1) - d(2)*sinTp)^2 < r(3)^2;
    cond2 = (x(1) - d(2)*cosTp)^2 + (y(1) + d(2)*sinTp)^2 > r(3)^2;
    if  ~(cond1 && cond2)
        x = NaN; y = NaN;
        return
    end

    %Step 5: Calculate the values of the coordinates of the relevant intersection points involving
    %circle 3
    xp13  =  (r(1)^2 - r(3)^2 + d(2)^2) / (2 * d(2));
    %yp13  = -0.5 / d(2) * sqrt( 2 * d(2)^2 * (r(2)^2 + r(3)^2) - (r(1)^2 - r(3)^2)^2 - d(2)^4 );
    yp13  = -0.5 / d(2) * sqrt( 2 * d(2)^2 * (r(1)^2 + r(3)^2) - (r(1)^2 - r(3)^2)^2 - d(2)^4 );

    x(2)   =  xp13*cosTp - yp13*sinTp;
    y(2)   =  xp13*sinTp + yp13*cosTp;

    xpp23 =  (r(2)^2 - r(3)^2 + d(3)^2) / (2 * d(3));
    ypp23 =  0.5 / d(3) * sqrt( 2 * d(3)^2 * (r(2)^2 + r(3)^2) - (r(2)^2 - r(3)^2)^2 - d(3)^4 );

    x(3) = xpp23*cosTpp - ypp23*sinTpp + d(1);
    y(3) = xpp23*sinTpp + ypp23*cosTpp;

end



function z = calcZoneAreas(nCircles, a, i)
    
    %Uses simple set addition and subtraction to calculate the zone areas
    %with circle areas a and intersection areas i

    if nCircles==2
        %a = [A1 A2]
        %i = I12
        %z = [A1-I12, A2-I12, I12]
        z = [a(1)-i, a(2)-i, i];
    elseif nCircles==3
        %a = [A1  A2  A3]
        %i = [I12 I13 I23 I123]
        %z = [A1-I12-I13+I123, A2-I12-I23+I123, A3-I13-I23+I123, ...
        %     I12-I123, I13-I123, I23-I123, I123];
        z = [a(1)-i(1)-i(2)+i(4), a(2)-i(1)-i(3)+i(4), a(3)-i(2)-i(3)+i(4), ...
                i(1)-i(4), i(2)-i(4), i(3)-i(4), i(4)];
    else
        error('')
        %This error gets caught earlier in the stack w. better error msgs
    end
end

function [Cx, Cy, aiz] = centroid2CI (x, y, r)

    %Finds the centroid of the area of intersection of two circles.
    %Vectorized to find centroids for multiple circle pairs
    %x, y, and r are nCirclePairs*2 arrays
    %Cx and Cy are nCirclePairs*1 vectors

    %Centroid of the area of intersection of two circles
    n = size(x,1);
    xic = zeros(n,2);
    az = zeros(n,2);
    
    dx = x(:,2)-x(:,1);
    dy = y(:,2)-y(:,1);
    d = sqrt(dx.^2 + dy.^2);
    
    %Translate the circles so the first is at (0,0) and the second is at (0,d)
    %By symmetry, all centroids are located on the x-axis.
    %The two circles intersect at (xp, yp) and (xp, -yp)
    xp = 0.5*(r(:,1).^2 - r(:,2).^2 + d.^2)./d;

    %Split the inner zone in two
    %Right side (Area enclosed by circle 1 and the line (xp,yp) (xp,-yp)
    %Angle (xp,yp) (X1,Y1) (xp,-yp)
    alpha = 2*acos(xp./r(:,1));
    %Area and centroid of the right side of the inner zone
    [xic(:,1) az(:,1)] = circleChordVals (r(:,1), alpha);
    %Angle (xp,yp) (X2,Y2) (xp,-yp)
    alpha = 2*acos((d-xp)./r(:,2));
    %Area and centroid of the left side of the inner zone
    [xic(:,2) az(:,2)] = circleChordVals (r(:,2), alpha);
    xic(:,2) = d - xic(:,2);

    %Thus the overall centroid  & area of the inner zone
    aiz = sum(az,2);
    Cx = sum(az.*xic,2)./aiz;
    
    %Now translate the centroid back based on the original positions of the
    %circles
    theta = atan2(dy, dx);
    Cy = Cx.*sin(theta) + y(:,1);
    Cx = Cx.*cos(theta) + x(:,1);
    
end

function centroidPos = zoneCentroids2 (d, r, Z)
    
    centroidPos = zeros(3,2);
    
    %Find the centroids of the three zones in a 2-circle venn diagram
    %By symmetry, all centroids are located on the x-axis.
    %First, find the x-location of the middle (intersection) zone centroid
    
    %Centroid of the inner zone
    centroidPos(3,1) = centroid2CI([0 d], [0 0], r);
    
    %Now, the centroid of the left-most zone is equal to the centroid of
    %the first circle (0,0) minus the centroid of the inner zone
    centroidPos(1,1) = -centroidPos(3,1)*Z(3)/Z(1);
    
    %Similarly for the right-most zone; the second circle has centroid at x=d
    centroidPos(2,1) = (d*(Z(2)+Z(3)) - centroidPos(3,1)*Z(3))/Z(2);
    
end

function centroidPos = zoneCentroids3 (x0, y0, d, r, Z)

    Z = Z(:);
        
    %Get area, points of intersection, and chord lengths
    [act, xi, yi, c, atr] = areaIntersect3Circ (r, d);
    atr = atr(:);
    r = r(:);
    
    %Area and centroid of the triangle within the circular triangle is
    xtr = sum(xi/3); 
    ytr = sum(yi/3);
        
    %Now find the centroids of the three segments surrounding the triangle
    i = [1 2; 1 3; 2 3]; 
    xi = xi(i); yi = yi(i);
    [xcs, ycs, acs] = circSegProps (r(:), x0(:), y0(:), xi, yi, c(:));
    
    %Overall centroid of the circular triangle
    xct = (xtr*atr + sum(xcs.*acs))/act;
    yct = (ytr*atr + sum(ycs.*acs))/act;
    
    %Now calculate the centroids of the three two-pair intersection zones
    %(Zones 12 13 23)
    %Entire zone centroid/areas

    %x, y, and r are nCirclePairs*2 arrays
    %Cx and Cy are nCirclePairs*1 vectors
    i = [1 2; 1 3; 2 3];
    [x2c, y2c, a2c] = centroid2CI (x0(i), y0(i), r(i));
    
    %Minus the three-circle intersection zone
    xZI2C = (x2c.*a2c - xct*act)./(a2c-act);
    yZI2C = (y2c.*a2c - yct*act)./(a2c-act);
    
    x0 = x0(:);
    y0 = y0(:);
    
    %Finally, the centroids of the three circles minus the intersection
    %areas
    i1 = [4 4 5]; i2 = [5 6 6];
    j1 = [1 1 2]; j2 = [2 3 3];
    x1C = (x0*pi.*r.^2 - xZI2C(j1).*Z(i1) - xZI2C(j2).*Z(i2) - xct*act)./Z(1:3);
    y1C = (y0*pi.*r.^2 - yZI2C(j1).*Z(i1) - yZI2C(j2).*Z(i2) - yct*act)./Z(1:3);
    
    %Combine and return
    centroidPos = [x1C y1C; xZI2C yZI2C; xct yct];
end


function [x, a] = circleChordVals (r, alpha)
    %For a circle centered at (0,0), with angle alpha from the x-axis to the 
    %intersection of the circle to a vertical chord, find the x-centroid and
    %area of the region enclosed between the chord and the edge of the circle
    %adapted from http://mathworld.wolfram.com/CircularSegment.html
    a = r.^2/2.*(alpha-sin(alpha));                         %Area
    x = 4.*r/3 .* sin(alpha/2).^3 ./ (alpha-sin(alpha));    %Centroid
end

function [xc, yc, area] = circSegProps (r, x0, y0, x, y, c)

    %Translate circle to (0,0)
    x = x-[x0 x0];
    y = y-[y0 y0];

    %Angle subtended by chord
    alpha = 2*asin(0.5*c./r);
       
    %adapted from http://mathworld.wolfram.com/CircularSegment.html
    area = r.^2/2.*(alpha-sin(alpha));                         %Area
    d   = 4.*r/3 .* sin(alpha/2).^3 ./ (alpha-sin(alpha));    %Centroid
   
    %Perpindicular bisector of the chord
    m = -(x(:,2)-x(:,1))./(y(:,2)-y(:,1));
    
    %angle of bisector
    theta = atan(m);
    
    %centroids
    xc = d.*cos(theta);
    yc = d.*sin(theta);
    
    %Make sure we're on the correct side
    %Point of intersection of the perp. bisector and the circle perimeter
    xb = (x(:,1)+x(:,2))/2;
    xc(xb<0) = xc(xb<0)*-1;
    yc(xb<0) = yc(xb<0)*-1;
    
    %Translate back
    xc = xc + x0;
    yc = yc + y0;
end


function [A0, I0, Z0, nCircles, fminOpts, vennOpts, patchOpts] = parseArgsIn (args)

    [A0, I0, Z0] = deal([]);
    nIn = length(args);
    badArgs = false;
    
    %Get the easy cases out of the way
    if nIn == 0
        badArgs = true;
    elseif nIn == 1
        %venn(Z)
        Z0 = args{1};
        nIn = 0;
    elseif nIn == 2
        if isnumeric(args{2})
            %venn (A,I)
            [A0, I0] = deal(args{1:2});
            nIn = 0;
        else
            %venn (Z, F)
            Z0 = args{1};
            args = args(2);
            nIn = 1;
        end
    else
        %Find the first non-numeric input arg
        i = find(~cellfun(@isnumeric, args), 1);
        if i == 2
            %venn(Z, ....)
            Z0 = args{1};
        elseif i == 3
            %venn(A, I, ...)
            [A0, I0] = deal(args{1:2});
        else
            badArgs = true;
        end
        nIn = nIn - i + 1;
        args = args(i:end);
    end
    
    if badArgs
        error('venn:parseInputArgs:unrecognizedSyntax', 'Unrecogized input syntax')
    end
    try
        [A0, I0, Z0] = parseInputAreas (A0, I0, Z0);
    catch
        error('venn:parseArgsIn:parseInputAreas', 'Incorrect size(s) for area vector(s)')
    end
    nCircles = length(A0);
    nZones = length(Z0);
             
    %Any arguments left?
    if nIn > 0 
        
        if isstruct(args{1})
            %FMIN search options
            f = args{1};
           
            nIn = nIn - 1;
            if nIn>0, args = args(2:end); end

            if length(f) == 1
                %Just double up
                fminOpts = [f f];
            elseif length(f) == 2
                %ok
                fminOpts = f;
            else
                error('venn:parseArgsIn', 'FMINOPTS must be a 1 or 2 element structure array.')
            end
        else
            %Use defaults
            fminOpts = [optimset('fminbnd'), optimset('fminsearch')];
        end
    else
        %Use defaults
        fminOpts = [optimset('fminbnd'), optimset('fminsearch')];
    end

    %If there's an even number of args in remaining
    if nIn>0 
        if mod(nIn, 2)==0
            %Parameter/Value pairs
            p = args(1:2:end);
            v = args(2:2:end);
            [vennOpts, patchOpts] = parsePVPairs (p, v, nZones);
        else
            error('venn:parseArgsIn', 'Parameter/Value options must come in pairs')
        end
    else
        vennOpts = defaultVennOptions;
        patchOpts = struct('Parameters', [], 'Values', []);
    end

end %parseArgsIn

function [vennOpts, patchOpts] = parsePVPairs (p, v, nZones)

    p = lower(p);

    %Break up P/V list into Venn parameters and patch parameters
    vennParamNames = {'plot', 'errminmode', 'parent'};
    [isVennParam, idx] = ismember(p, vennParamNames);
    idx = idx(isVennParam);
    %vennParams = p(isVennParam);
    vennVals = v(isVennParam);
    
    %First do Patch options
    patchOpts.Parameters = p(~isVennParam);
    patchOpts.Values = v(~isVennParam);
        
    %Now do Venn options
    vennOpts = defaultVennOptions;
        
    %PLOT
    i = find(idx==1, 1);
    if i
        plot = lower(vennVals{i});
        if islogical(plot)
            vennOpts.Plot = plot;
        else
            if ischar(plot) && any(strcmp(plot, {'on', 'off'}))
                vennOpts.Plot = strcmp(plot, 'on');
            else
                error('venn:parsePVPairs', 'Plot must be ''on'', ''off'', or a logical value.')
            end
        end
    end
    
    %ERRMINMODE
    i = find(idx==2, 1);
    if i
        mode = lower(vennVals{i});
        okModes = {'None', 'TotalError', 'ChowRodgers'};
        [isOkMode, modeIdx] = ismember(mode, lower(okModes));
        if isOkMode                
            vennOpts.ErrMinMode = okModes{modeIdx};
        else
            error('venn:parsePVPairs', 'ErrMinMode must be None, TotalError, or ChowRodgers')
        end
    end

    %PARENT
    i = find(idx==5, 1);
    if i
        h = v{i};
        if length(h)==1 && ishandle(h) 
            vennOpts.Parent = h;
        else
            error('venn:parsePVPairs', 'Parent must be a valid scalar handle')
        end
    end
    
       
end %parsePVPairs

function [A0, I0, Z0] = parseInputAreas (A0, I0, Z0)

    %Switch to row vectors
    A0 = A0(:)';
    I0 = I0(:)';
    Z0 = Z0(:)';

    if isempty(Z0)
        %A0 and I0 supplied
        
        Z0 = calcZoneAreas (length(A0), A0, I0);
    else
        %Z0 supplied
        switch length(Z0)
            case 3
                A0 = Z0(1:2)+Z0(3);
                I0 = Z0(3);
            case 7
                A0 = Z0(1:3)+Z0([4 4 5])+Z0([5 6 6])+Z0(7);
                I0 = [Z0(4:6)+Z0(7) Z0(7)];
            otherwise
                error('')
        end
    end
end

function vennOpts = defaultVennOptions 
    
    vennOpts = struct(...
        'Plot'          ,true               ,...
        'Labels'        ,[]                 ,...
        'PopLabels'     ,false              ,...
        'DrawLabels'    ,false              ,...
        'Parent'        ,[]                 ,...
        'Offset'        ,[0 0]              ,...
        'ErrMinMode'    ,'TotalError'       );
    
end

function [d, x, y, A, I, Z] = preallocVectors (nCirc)

    %Initialize position vectors
    x = zeros(1, nCirc);
    y = zeros(1, nCirc);

    if nCirc==2
        d = 0;
        I = 0;    
        A = zeros(1,2);
        Z = zeros(1,3);

    else %nCirc==3
        d = zeros(1,3);
        I = zeros(1,4);
        A = zeros(1,3);
        Z = zeros(1,7);
    end
end

    
%% Another internal function 

function colourBoxPlot(plotData,groups, color, includeScatter)

% set the color to the box plots
if nargin == 2 || isempty(color)
    rng(6);
    color = rand(length(unique(groups )),3) ;
end

% plot the data
figure()
boxplot(plotData,groups,'Color', flipud(color) ,'Symbol','k+', ...
    'OutlierSize',5) ;

% set some figure properties and add title ot the figure
set(gca,'FontSize',14,'LineWidth',1.5,'Box','off')

% set the line width of the box plots
set(findobj(gca,'type','line'),'linew',2)
set(findobj(gca,'Tag','Lower Whisker'),'LineStyle','-')
set(findobj(gca,'Tag','Upper Whisker'),'LineStyle','-')

% set the color of the box plots
h4 = findobj(gca,'Tag','Box') ;
for kk=1:length(h4)
    patch(get(h4(kk),'XData'),get(h4(kk),'YData'),...
        color(kk,:),'FaceAlpha', 0.3,'LineStyle','-');
end

% add a scatter plot if we that is true
if includeScatter
   
    % get the unique groups 
    uniqueGroups = unique(groups) ;

    % add the scatter plots
    hold on
    
    % add scatter plot to the box plots
    groupSc1 = plotData(groups == uniqueGroups(1,1));
    groupSc2 = plotData(groups == uniqueGroups(2,1));
    
    % set the size the size of the marker
    makerSz = 30;
     
    % get only a smaller subset if the data points are too many
    if length(groupSc1) > 600
        groupSc1 = randsample(groupSc1,600) ;
        
        % change the markerSize 
        makerSz = 15 ;
    end
    
    if length(groupSc2) > 600
        groupSc2 = randsample(groupSc2,600) ;
    end
    
    x = ones(length(groupSc1)).*(1+(rand(length(groupSc1))-0.5)/5) ;
    x1 = ones(length(groupSc2)).*(1+(rand(length(groupSc2))-0.5)/10);
    
    % here is the first scatter plot
    scatter(x(:,1), groupSc1, makerSz, color(2,:),'filled', ...
        'MarkerFaceColor',color(2,:),'Marker','o','MarkerFaceAlpha',0.8)
    hold on
    scatter(x1(:,2).*2, groupSc2, makerSz, color(1,:),'filled', ...
        'MarkerFaceColor',color(1,:),'Marker','o','MarkerFaceAlpha',0.8)
    
    hold off
end

end

%% Another internal function 

function ManhattanPlot( filename, varargin )

%   This function takes a GWAS output file  and plots a Manhattan Plot.

%%  ARGUMENTS:
%   sex: defaults to 0, set to 1 to include sex chromosomes

%   sig: defaults to 5e-8, significance threshold for horizontal line

%   vert: defaults to 0, set to 1 for tick labels have chr in them and are
%   written upwards

%   labels: defaults to [-1,-1]. Set to [x,y] to label the top SNP on each
%   locus with p<x. Locus defined by windows of y base pairs.

%   outfile: defaults to filename of input file, set to something else to
%   change name of output file

%   title: defaults to 'Manhattan Plot'. Fairly self-explanatory

%   save: defaults to 1, set to 0 to disable saving to save time

%   format: defaults to PLINK, options {PLINK,BOLT-LMM,SAIGE} This is used
%   to identify the correct column headings. If using anything else, rename
%   the header line so that CHR, BP, and P are the headers for chromosome,
%   base pair and p value, and the default PLINK should catch it

%%   Usage:
%   ManhattanPlot('gwas.assoc.fisher',varargin) will take the association
%   analysis in 'gwas.assoc.fisher, generate a Manhattan Plot, and store it
%   in gwas_ManhattanPlot.png, which should be publication-ready, and
%   gwas_ManhattanPlot.fig for minor readjustments in MATLAB's GUI. This is
%   fine as .fisher is a PLINK format file. For a SAIGE output (and
%   significance threshold 5e-9) with sex chromosomes shown, and the top
%   hit for each SNP within 1 mb and below p=1e-6 labelled, use
%   ManhattanPlot('gwas.saige','format','SAIGE','sig',5e-9,'sex',1,
%   'labels',[1e-6,1000000])

%   Tested using assoc, assoc.logistic and assoc.fisher files generated by
%   Plink 1.7. Also tested using BOLT-LMM and SAIGE.
%
%   Harry Green (2019) Genetics of Complex Traits, University of Exeter

%% Reading in optional arguments and input file

p = inputParser;
defaultsex = 0;
defaultvert = 0;
defaultsig = 5e-8;
defaultsave = 1;
defaultoutfile = filename;
defaulttitle= 'Manhattan Plot';
defaultlabels = [-1, -1];
defaultformat = 'PLINK';
expectedformat = {'PLINK','BOLT-LMM','SAIGE'};

addRequired(p,'filename');
% addOptional(p,'outfile',defaultoutfile,@ischar);
addOptional(p,'sex',defaultsex,@isnumeric);
addOptional(p,'vert',defaultvert,@isnumeric);
addOptional(p,'labels',defaultlabels,@isnumeric);
addOptional(p,'sig',defaultsig,@isnumeric);
addOptional(p,'save',defaultsave,@isnumeric);
addParameter(p,'format',defaultformat,...
    @(x) any(validatestring(x,expectedformat)));
addParameter(p,'outfile',defaultoutfile,@ischar);
addParameter(p,'title',defaulttitle,@ischar);

parse(p,filename,varargin{:});

filename=p.Results.filename;
sex=p.Results.sex;
vert=p.Results.vert;
sig=p.Results.sig;
save=p.Results.save;
format=p.Results.format;
outfile=p.Results.outfile;
labels=p.Results.labels;
plottitle=p.Results.title;

% ################ Check if the input is a table for a file name #########

if ischar(filename)
    % MATLAB refuses to read file formats like .assoc, so first rename as
    % .txt
    try
        copyfile(filename, strcat(filename,'.txt'));
        opts = detectImportOptions(strcat(filename,'.txt'),...
            'NumHeaderLines',0);
        % the columns of T will be the headers of the assoc file
        T = readtable(strcat(filename,'.txt'),opts);
        
        %delete unwanted .txt file that we didn't want anyway
        delete(strcat(filename,'.txt'))
    catch
        T = readtable(filename);
    end
elseif istable(filename)
    T = filename ;
end

%% File format defintitions
% This section uses the default column headings from PLINK, BOLT-LMM and
% SAIGE. Easy to modify for other software packages

% check the CHR variable exist then change to PLINK format
if ~ismember(T.Properties.VariableNames,'CHR')
    T.Properties.VariableNames([1,2]) = {'CHR','BP'} ;
end

if strcmp(format,'PLINK')
    % These are the only useful columns for a manhattan plot, chromosome,
    % base pair, p value
    try 
        tab2=[T.CHR,T.BP,T.P]; 
    catch
        % some time this does not work
        locBP = find(ismember(T.Properties.VariableNames,'POS'), true);
        T.Properties.VariableNames(locBP) = "BP";
        tab2=[T.CHR,T.BP,T.P]; 
    end
end

if strcmp(format,'BOLT-LMM')
    % These are the only useful columns for a manhattan plot, chromosome,
    % base pair, p value
    tab2=[T.CHR,T.BP,T.P_BOLT_LMM]; 
end

if strcmp(format,'SAIGE')
    % These are the only useful columns for a manhattan plot, chromosome,
    % base pair, p value
    tab2=[T.CHR,T.POS,T.p_value]; 
end

%%

if sex==0
    tab2=tab2(tab2(:,1)<23,:); %remove chr23 if not wanting sex chromosomes
end

% variable to track base pairs, this helps the plotting function to know
% where to start plotting the next chromosome
bptrack=0; 
tab2(tab2(:,3)==0,:)=[];
for i=1:22+sex
    hold on
    %a scatterplot. On the x axis, the base pair number + bptrack, which
    %starts at 0. On the y axis, - log 10 p
    plot(tab2(tab2(:,1)==i,2)+max(bptrack),-log10(tab2(tab2(:,1)==i,3)),'.'); 
    
    %this updates bptrack by adding the highest base pair number of the
    %chromosome. At the end, we should get the total number of base pairs
    %in the human genome. All values of bptrack are stored in a vector.
    %They're useful later
    bptrack=[bptrack,max(bptrack)+max(max(tab2(tab2(:,1)==i,:)))]; 
end

%if strongest hit is 1e-60, plot window goes up to 1e-61
ylimit=max(max(-log10(tab2(:,3)))+1); 
ylimitmin=min(min(-log10(tab2(:,3)))); %and down to the highest p value.

%genome wide significant line, uses the sig optional argument which
%defaults to 5e-8
plot([0,max(bptrack)],-log10([sig,sig]),'k--') 

% ylabel('$-\log_{10}(p)$','Interpreter','latex')
ylabel('-log_1_0(p)')
xlim([0,max(bptrack)])
%y axis will always go to the significance threshold/5 at least
ylim([floor(ylimitmin),max(ylimit,-log10(sig/5))])

%this calculates the moving average of bptrack and uses them as chromosome
%label markers. This puts them in the middle of each chromosome.
M=movmean(bptrack,2); 
xticks(M(2:end));

%% Rotation section
% this section of the code changes the x axis label orientation.

if vert ==0
    xlabel('Chromosome')% ,'Interpreter','latex')
    xticklabels( 1:23 );
end
if vert ==1
    xtickangle(90)
    xticklabels( {'chr1','chr2','chr3','chr4','chr5','chr6','chr7',...
        'chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr15',...
        'chr16','chr17','chr18','chr19','chr20','chr21','chr22','chrXY'});
end

%% Annotation section
% this section of the code annotates the top SNP on each locus

% Neither of these should be negative. Defaults to -1 -1
if ~(labels(1)==-1||labels(2)==-1) 
    labellim=labels(1);
    labelprox=labels(2);
    
    %as nothing else will be used for labelling, the reduced table is
    %easier to manage
    tab3=tab2(tab2(:,3)<labellim,:); 
    
    %easier if they're in ascending order of p value;
    [~,index]=sort(tab3(:,3));
    tab3=tab3(index,:);
    
    %this becomes a list of SNPs that have been labelled. It starts with a
    %dummy SNP to give the first SNP something to compare with. This will
    %always fail so the first label always gets added
    labelledsnps=[0,0,0]; 
    
    for i=1:size(tab3,1)
        test=abs(tab3(i,:)-labelledsnps);
        %if there are no snps on the same chromosome and within 1 MB
        if sum(test(:,1)==0&test(:,2)<labelprox)==0 
            
            %this plots a square over the corresponding dot
            plot(bptrack(tab3(i,1))+tab3(i,2),-log10(tab3(i,3)),'ks',...
                'MarkerSize',10,'LineWidth',1) 
            
            %this puts together a strong of chrx:y
            labels=strcat('chr',string(tab3(i,1)),':',...
                string(bptrack(tab3(i,1))+tab3(i,2))); 
            
            %this plots the label above and to the left of the dot, shifted
            %up by 0.05 to give some space
            text(bptrack(tab3(i,1))+tab3(i,2),-log10(tab3(i,3))+0.05,...
                char(labels),'VerticalAlignment','bottom',...
                'HorizontalAlignment','left') %,'Interpreter','latex') 
            labelledsnps=[labelledsnps;tab3(i,:)];
        end
    end
end
%% finalising file

% takes the output file name up to the first decimal point
a = strsplit(outfile,'.'); 

% my code 
set(gca,'FontSize',12,'LineWidth',1)

box on
title(plottitle,'fontsize',14)
% title(plottitle,'Interpreter','latex','fontsize',14)

if save~=0
    name=[a{1},'_ManhattanPlot.fig']; %adds _ManhattanPlot to filenames
    savefig(name);
    set(gcf, 'paperunits', 'centimeters', 'paperposition', [0 0 30 10])
    print([a{1},'_ManhattanPlot.png'],'-dpng','-r300', '-painters') % -r600
end



end
        
% ******************* Another Internal function ******************

function exitPos = ultraBars(inData, colors, rowNames, ...
    horizonBars , nextAxis)

% Input:
% inData: a matrix and vector of plotting data
% colors: colors for each unique value of inData
% rowNames: a cell array for name for each row
% horizonBar: specifies whether to add a horizontal bar to each plot or not
 
% get the number of unique numbers and remove the 0 which means not plot
uniqueVars = unique(inData) ;
uniqueVars(uniqueVars == 0) = [] ;

% create the legend variable if the inData is a row vector
if size(inData,1) == 1
    lgdVar = split( num2str(uniqueVars) )';
end

% get the number of colours to plot
if ~exist('colors','var') || isempty(colors)
    colors = rand(length(uniqueVars), 3);
end

% check the are is a row name for each row in inData
if size(inData,1) > 1 % only valid for matrix data
    if size(inData,1) ~= size(rowNames,1)
        error('row names should contain a name for each row in plotting data')
    end
end

% check if the orizontal bars are required on the plot
if ~exist('horizonBars','var')
    horizonBars = false;
end

% check for color errors 
if length(uniqueVars) > size(colors,1) 
    error('A color must be specified for each value in Data')
end

% plot the figure for row vectors 
% I have separated the two peaces of code because one does not need a box
% around the data
if size(inData,1) == 1
    % figure()
    % axes('position',[0.1300,0.11,0.7,0.04]);
    % set the bar width for large case
    if size(inData,2) < 500
        barwidth = 0.9 ;
    else
        barwidth = 1 ;
    end
    
    % now plot the data
    for ii = 1:length(uniqueVars)
        % get the data for that values and plot it in that color
        plotData = double(ismember(inData,uniqueVars(ii) )) ;
        bar(plotData,'FaceColor',colors(ii,:),...
            'EdgeColor',[1 1 1] ,'BarWidth',barwidth) ;
        hold on
    end
    set(gca,'GridColor',[1 1 1], 'XLim', [0.5 size(inData,2)+0.5 ], ...
        'XColor',[1 1 1] ,'YColor',[1 1 1],'XTick',[] , 'YTick',[],...
        'FontWeight','bold')
    % add a legend to the bar using the legendflex function
   
% if the is more thhan one row in the dataset
else
    % make a plot of multiple bar chart of top of each other
    % add the subtype plots to the heatmap using a loop
    % initialise some variables
    global plotTime
    figure(plotTime*100)
    % set(gcf,'position',[100,50,800,600])
        
    yInitial = 0.05; yPosSaved = yInitial;
    ySize = 0.02; increaseby = 0.022; % 0.44 0.4
    % change the size the bar if plotTime == 3
    if plotTime == 3
        ySize = 0.05; increaseby = 0.055;
    end
    xEndPos = 0.7 ; % 0.7750
    % loop over the rows and ascend by column
    for jj = 1:size(inData,1) % begin plotting from the bottem
        % define the next axis for the follow up plots
        if ~exist('nextAxis','var') || isempty(nextAxis)
            axes('position',[0.1300,yInitial,xEndPos,ySize]);
        elseif exist('nextAxis','var') && jj > 1
            axes('position',[0.1300,yInitial,xEndPos,ySize]);   
        else
            axes('position',nextAxis);
            yInitial = nextAxis(2) ;
            ySize = nextAxis(4) ; xEndPos = nextAxis(3) ;
            yPosSaved = yInitial ;
        end         
        for ii = 1:numel(uniqueVars) 
            plotData = double(ismember(inData(jj,:),uniqueVars(ii) )) ;
            bar(plotData,'FaceColor',colors(ii,:),'EdgeColor',[1 1 1] ,...
                'BarWidth',0.9) ;   
            hold on
            
            % add the name of the genes to the left of heatmap
            if exist('rowNames','var')
                dim = [0.02 yInitial 0.11 increaseby];
                annotation('textbox',dim,'String',rowNames{jj},...
                'FitBoxToText','on','FontSize',10,'EdgeColor','none',...
                    'HorizontalAlignment','right','FontWeight','bold',...
                    'VerticalAlignment','middle');
            end
        end
        % change the plot properties
        set(gca,'GridColor',[1 1 1], 'XLim', [0.5  size(inData,2)+0.5],...
            'XColor',[1 1 1] ,'YColor',[1 1 1],'YTickLabel',[],...
            'XTickLabel',[],'FontWeight','bold','YTick',[],'XTick',[])
        % increase the value to change the colors and plot positions
        yInitial = yInitial + increaseby;
    end

    % add a grey box to the plot usign the annotation function
    if plotTime ~= 3 % dont add the box if plotTime is == 3
        dim = [0.1300, yPosSaved, xEndPos, increaseby*size(inData,1)];
        annotation('rectangle',dim ,'Color',[0.5, 0.5, 0.5])
    end
    hold off 
   
    % plot the horizontal bars if they are required
    % prellocate the bar data size
    barhData = zeros(size(inData,1),numel(uniqueVars)) ;
    if horizonBars == true
        axes('position',[xEndPos+0.137, yPosSaved, 0.12, ...
            increaseby*size(inData,1)+0.001 ]);
        for kk = 1:numel(uniqueVars) 
            barhData(:,kk) = sum(inData == uniqueVars(kk),2) ;   
        end
        bar1 = barh(barhData,'stacked','BarWidth',0.85) ;
        % make sure there are no colors and spaces between the axis and the
        % first and last bar
        set(gca,'GridColor',[1 1 1], 'YLim', [0.5 size(inData,1)+0.52 ], ...
            'XColor',[0.5 0.5 0.5] ,'YColor',[1 1 1],'FontSize',10,...
            'YTick',[],'FontWeight','bold','Box','off','TickDir', 'out',...
            'LineWidth',1,'XAxisLocation','origin')
        
        % annoate the bar graph
        for ii = 1:size(barhData,2)
            set(bar1(ii),'FaceColor',colors(ii,:))
        end
        
    end
    exitPos = [0.1300,yInitial,xEndPos, ySize] ;
end

end

