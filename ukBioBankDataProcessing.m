% Process UK BioBank Clinical Data
% This script returns only the white and black data from the UK BioBank
% data 

% fresh start
clc; close; clear 

% change directory
cd('/Users/sinkala/Documents/MATLAB/UKBioBank')

%% Load the dictionary data and Show Case data 

% first the dictionary
theDict = readtable('Data_Dictionary_Showcase.xlsx') ;

% how many samples have complete data
summary(categorical(theDict.Stability))

% return only data for which the collection was complete by the time that
% we got these data data from the UK biobank
theDict = theDict( ismember(theDict.Stability,'Complete'), :);

% load the showcase data 
showCase = readtable('Codings_Showcase.xlsx');

% add the "FieldId" from the theDict to the show case data 
showCase = innerjoin( theDict(:,{'FieldID','Field','Coding'}), showCase);

% lets see the ethic groups
ethnicGroups = showCase(showCase.FieldID == 21000, :) ;
disp(ethnicGroups)

%% Load the all the UK Biobank data 

% create a datastore
ds = datastore('ukb40800.csv');

% Preview the data

% After creating the datastore, you can preview the data without having to
% load it all into memory. You can specify variables (columns) of interest
% using the SelectedVariableNames property to preview or read only those
% variables.

% % get the variable names related to the ethic group (1002)
selectedGrp = ds.SelectedVariableNames(contains( ...
    ds.SelectedVariableNames,{'eid','21000','21001','21002','21003'}));

% lets preview those data
ds.SelectedVariableNames = selectedGrp  ;
preview(ds)

%% ====================== NOTE about the data =========================

% ** It turn out that the FieldID in the Data_Dictionary_Showcase.csv
% files is what is used in the ukb40800.csv as the variable names.

% In the  Data_Dictionary_Showcase.csv the stability column can be used to
% select variable that are "Complete". I can leave out those that are still
% "Accruing"

% In the  Data_Dictionary_Showcase.csv the show the number of instance the
% particular variablename appear in the ukb40800.csv data. i.e. 4 mean that
% a particular variablename appear 4 times in the ukb40800.csv data. e.g.,
% the "21000" shows 3 instance and it appear in the ukb4000.csv as
% 21000_0_0    21000_1_0    21000_2_0. However, I seem that the 0_0
% variable is the one with complete data

% Some of the FieldID(s) (variablenames) from the
% Data_Dictionary_Showcase.csv file do not appear in the ukbiobank.csv
% files. For example 20160 does not appear.

% In the Data_Dictionary_Showcase.csv file the variable name coding is
% matched to the Codings_Showcase.csv files. e.g., the coding for "21000"
% is "1001" is values is avaiable in both files.

% The "Value" variable is the Codings_Showcase.csv file indicates variable
% under a specific "Coding" related to a particular "FieldID". e.g., Under
% the "FieldID" "21000" which represents ethic group, the "coding" of
% "1001" in the has many rows that mean different ethic groups in the
% Codings_Showcase.csv

%% Process the data

% reload the datastore
ds = datastore('ukb40800.csv');

% convert everything to a string so that I can be read in MATLAB -- This
% was the major problem
ds.SelectedFormats = repmat({'%q'}, 1, length(ds.SelectedFormats) );

% only select variable that have a 0_0
% get all the variable names from the dictionary and include the sample ID
% which I missed out in the first iteration 
dictVars = [ds.VariableNames(1) ; ...
    strcat('x',string(theDict.FieldID) ,'_0_0') ];
ds.SelectedVariableNames = ds.SelectedVariableNames( ismember( ...
            ds.SelectedVariableNames, dictVars ));
        
% preallocate the data table
allData = [];  

% Get a number of partitions to parallelize datastore access over the
% current parallel pool.
n = numpartitions(ds,gcp);

% time the running of this programme
tstart = tic ;

parfor kk = 1:n
    
    % Partition the datastore and read the data in each part.
    subds = partition(ds,n,kk);
    
    % By calling the read function within a while loop, you can perform
    % intermediate calculations on each subset of data, and then aggregate
    % the intermediate results at the end.
    while hasdata(subds)
        
        % read in a smaller dataset: note that this return only the data
        % that has complete record and those that were taken at the first
        % instance i.e., have _0_0
        curData = read(subds);
        
        % add the sample ID variable names
        curData.Properties.VariableNames(1) = "SampleID" ;
        
        % clear up the variable names
        curData.Properties.VariableNames(2:end) = extractBetween( ...
            curData.Properties.VariableNames(2:end), 'x','_') ;
        
        % insert the codes that are found in the ShowCase files into the
        % curData table
        for ii = 2:width(curData)
            
            % get the ShowCase data that matched the current data column
            % (variable)
            curCases = showCase(showCase.FieldID == ...
                str2double(curData.Properties.VariableNames(ii)), :) ;
            
            % check if the current Field has a code in the showCase file
            if isempty(curCases)
                continue
            end
            
            % if the "Value" variable is a double in the showCase file,
            % then covert it to a string
            if isreal(curCases.Value)
                curCases.Value = string(curCases.Value) ;
            end
            
            % change the coding in the table to that present in the
            % showCase file
            [these, locCodes] = ismember(curData.(ii),curCases.Value);
            
            % remove the zero values
            locCodes(locCodes == 0) = [] ;
            
            % find where non-zero values these order in the curData
            these = find(these);
            
            % now add the values to the table
            curData.(ii)(these) = curCases.Meaning(locCodes) ;
            
        end
        
        % change the codes contained in each column of the data to codes in
        % instances given in the showCase files.
        
        % add the correct variable names and covert the data to the correct
        % format based on the what the Data_Dictionary_Showcase.csv file
        % says
        for ii = 2:width(curData)
            
            % get the current variable names
            curVarName = str2double(curData.Properties.VariableNames(ii));
            
            % replace the current variable name which is a number with the
            % name from the theDict.Field
            newVarName =  matlab.lang.makeValidName( ...
                theDict.Field(theDict.FieldID == curVarName) ) ;
            
            % check that the variable names exists in the curData table
            if isempty(newVarName)
                continue
            end
            
            % check if the variable names already exist in the table
            if any(ismember( curData.Properties.VariableNames(1:ii), ...
                    newVarName ) )
                % find the  number of instances that the variable names
                % present in the table and use that make a number variable
                % name
                numOfOccurs = sum( ismember( ...
                    curData.Properties.VariableNames(1:ii),newVarName )) ;
                
                % create a new string
                newVarName = strcat(newVarName, string(numOfOccurs)) ;
            end
            
            % then change the varible name
            curData.Properties.VariableNames(ii) = newVarName;
            
            % ====== convert the variable to correct format ====
            
            % check the types format this variable should be
            ValueType = theDict.ValueType(theDict.FieldID == curVarName) ;
            
            if isempty(ValueType)
                continue
            end
            
            % now do the conversion
            if contains(ValueType,'Categorial','IgnoreCase',true)
                % convert to a categorical array
                curData.(ii) = categorical( curData.(ii)) ;
            elseif ismember(ValueType,'Continuous') || ...
                    ismember(ValueType,'Integer')
                % convert to a double
                curData.(ii) = str2double( curData.(ii)) ;
            end
        end
        
        % concatenate the data to a table
        allData = [allData ; curData ] ;
        
        telapsed = toc(tstart)/60 ;
        
        fprintf('\n Time elapsed is %s Minutes \n', num2str(telapsed) )
        
    end
end

% return only data of blacks and white and save the data of only black and
% whites
allDataBnW = allData( ismember( allData.EthnicBackground, ...
    {'African','British','White','Irish','Any other white background'}),:);

% add a new variable name that combineds the British and the Whites into
% one group
myEthicGroup = mergecats( categorical(allDataBnW.EthnicBackground), ...
    {'White','British','Irish','Any other white background'},'White') ;
allDataBnW = addvars( allDataBnW, myEthicGroup, 'NewVariableNames', ...
    'EthnicGroupBW','Before',1) ;

% ===================== Clean up the variable names ================
% remove participants that have withdrawn from the UK Biobank study 

% load the Ids of participants that have withdrawn from the study
ids = readtable('w53163_20210201.csv'); % from Nicky
ids2 = readtable('w53163_20200820.csv'); % from Samar

% put the ids together 
ids = [ids; ids2] ;
ids.Var1 = cellstr(num2str(ids.Var1));

% get only the data of the participants that remain in the study and save
% to excel 
allDataBnW = allDataBnW(~ismember( allDataBnW.SampleID, ids.Var1), :) ;
fprintf('\n Now saving the BnW data to .csv file\n')
tic
writetable(allDataBnW,'cleaned UK Biobank Blacks and Whites.csv')
toc

% get only the data of the participants that remain in the study and save
% to excel 
allData = allData(~ismember(allData.SampleID, ids.Var1), :) ;
fprintf('\n Now saving the all the data to .csv file\n')
tic
writetable(allData,'all_cleaned_UK_Biobank.csv')
toc

writetable('hpc_ukbiobank_clinical.csv')

% remove the long name
data = allDataBnW ;

% % lets see the data 
% data.SampleID = str2double(data.SampleID ) ;
% data = data(:,[1,2]);

% remove some of the the variables
clear ii ValueType curData numOfOccurs newVarName curVarName these ...
    locCodes numOfTime selectedGrp ans curCases dataSegs fileExt ...
    numOfTimes  myEthicGroup allData dictVars timIDs allDataBnW ...
    D a ans barbase freq hBar dataAll ds ids ids2 n startTime telapsed ...
    tstart

%% save the data in smaller chucks 

% if the dataset is too large -- which it is not 
    
% % create an array for the chucks and the number to append to the filename 
% dataSegs = 1:round(height(allDataBnW)/10):height(allDataBnW) ;
% 
% % create a director for the cleaned data 
% if ~exist('cleanedData','Dir')
%     mkdir cleanedData
%     addpath /Users/sinkala/Documents/MATLAB/UKBioBank/cleanedData
%     cd /Users/sinkala/Documents/MATLAB/UKBioBank/cleanedData
% else
%     cd /Users/sinkala/Documents/MATLAB/UKBioBank/cleanedData
% end
% 
% % save the data by segments
% for ii = 1:length(dataSegs)-1
%    
%     % print something to the screen 
%     fprintf(['\n Saving data to .csv file # %d for the range ', ...
%          'between %d to %d \n'],ii, dataSegs(ii), dataSegs(ii+1) )
%     
%     % get a smaller segment of the daata without adding the previous row to
%     % the data 
%     if ii == 1 
%         curData = allDataBnW(dataSegs(ii):dataSegs(ii+1),:) ;
%     else
%         curData = allDataBnW(dataSegs(ii)+1:dataSegs(ii+1),:) ;
%     end
%     
%     tic
%     
%     % save that segment to a csv file
%     writetable(curData,sprintf('cleaned UK Biobank_%d.csv',ii));
%     
%     toc
%     
% end
% 
% clear ii ValueType curData numOfOccurs newVarName curVarName these ...
%     locCodes numOfTime selectedGrp ans curCases dataSegs fileExt ...
%     numOfTimes 

%% Run some T-test and chi-Square test to find interesting variables

% these variable should look at conditions that show a signficant
% difference between blacks and white. Based of the p-values of the t-test
% for continous variables and chi-square test for other types of variables.

% specify the directory where the results will be save on the computer and
% create that folder on the computer
myDir = '/Users/sinkala/Documents/MATLAB/UKBioBank/projectFigures' ;
if ~exist(myDir,'dir')
    mkdir(myDir)
end

% % load the data
if ~exist('data','var')
    data = readtable('cleaned_UK_Biobank_Blacks_and_Whites.csv');
end

% change the variable name 
data.Properties.VariableNames{1} = 'EthnicGroupBW';

if iscell(data.EthnicGroupBW)
   data.EthnicGroupBW = categorical(data.EthnicGroupBW);
end

% preallocate the table for the t-test and chi-square test results
ttestTable = table(cell(0,0));
chiTable = table(cell(0,0));

% run the data in a loop
for ii = 3:width(data)
    
    % check that the data in the column are intergers then perform a t-test
    if isnumeric(data.(ii))
        
        % print something to the screen
        if rem(ii,50) == 0
            fprintf('\nRunning t-test # %d of %d Variables \n', ...
                ii, width(data)-2)
        end
        
        % add the variable name to the t-test table
        ttestTable(ii,1) = data.Properties.VariableNames(ii) ;
        
        % get the data in the current column and fill the outliers using
        % linear interpolation 
        curCol = filloutliers(data.(ii),'linear') ;
        
        % perform a ttest with unequal variable assumed
        [~,p,ci,stats] = ttest2(curCol(data.EthnicGroupBW =='African'),...
            curCol(data.EthnicGroupBW == 'White'),'Vartype','unequal') ;
        
        % get the mean values for the two groups
        meanAfricans = nanmean(curCol(data.EthnicGroupBW == 'African')) ;
        meanWhites = nanmean(curCol(data.EthnicGroupBW == 'White')) ;
        
        % add to the table
        ttestTable(ii,2:7) = num2cell([meanAfricans,meanWhites, ...
            stats.tstat,ci(1), ci(2) , p]);
        
        % also plot a figure and save the figure
        plotName = replace(data.Properties.VariableNames{ii}, '_','-');
        xValues = curCol ;
        groups = data.EthnicGroupBW=='African';
        createScatterBoxPlot( xValues, groups, p , plotName) 
        
        % save the figure
        saveas(gcf,[myDir,'/',plotName,'.png'],'png');
        
        % close the figure
        close ;
        
    else % are data either categorical or cell arrays (tuples)
        % print something to the screen
        if rem(ii,50) == 0
            fprintf('\nRunning chi-square test # %d of %d Variables \n',...
                ii, width(data)-2)
        end
          
        % peform the chi-square test
        if iscategorical(data.(ii)) % for categorical data 
            % get the data and remove the rows with missing data 
            curData = addvars( data(:,1), data.(ii), ...
                'NewVariableNames','newVar') ;
            curData( ismissing(curData.newVar), :) = [];
            
            % perform the chi-square test
            [~,chi2,chiPvalue] = crosstab(curData.EthnicGroupBW,...
                curData.newVar) ;
            
        else % for cell arrays (tuples)
            % try to convert the variables to a categorical array and if it
            % fails then continue
            curData = addvars( data(:,1), categorical(data.(ii)), ...
                'NewVariableNames','newVar') ;
            % remove the rows with missing categoroes
            curData(ismissing(curData.newVar), :) = [];
            
            % now perform the chi-square test
            try
                [~,chi2,chiPvalue] = crosstab(...
                    curData.EthnicGroupBW, curData.newVar ) ;
            catch
                continue
            end
        end
        
        % add the variable name to the table
        chiTable(ii,1) = data.Properties.VariableNames(ii) ;
        
        % add the statistics to the table 
        chiTable(ii,2:3) = num2cell([chi2,chiPvalue]);
    end
end
    
% add the variable names to the table and then remove the rows with missing
% results and add the adjusted p value to the table 
ttestTable.Properties.VariableNames = ...
    {'Variable','meanAfricans','meanWhites','tValue','lowerBound',...
    'upperBound','pValue'} ;
ttestTable( cellfun(@isempty, ttestTable.Variable ), :) = [] ;
ttestResults = addvars( ttestTable, ...
    mafdr(ttestTable.pValue,'BHFDR',true),'After','pValue', ...
    'NewVariableNames','FDR') ;
ttestResults = sortrows(ttestResults,'pValue','ascend');

% add the variable names to the table and then remove the rows with missing
% results and add the adjusted p value to the table 
chiTable.Properties.VariableNames = {'Variable','chiValue','pValue'};
chiTable( cellfun(@isempty, chiTable.Variable ), :) = [] ;
chiResults = addvars( chiTable, ...
    mafdr(chiTable.pValue,'BHFDR',true),'After','pValue', ...
    'NewVariableNames','FDR') ;
chiResults = sortrows(chiResults,'pValue','ascend');

% save the results to excel 
writetable(chiResults,'Group Comparisons Results.xlsx','Sheet','ChiSquare')
writetable(ttestResults,'Group Comparisons Results.xlsx','Sheet','T-test')

clear ii curData tbl chi2 chiPvalue labels ci p stats meanAfricans  ...
    meanWhites chiTable ttestTable chiData xValues groups plotName ...
    curCol

%% ================== The Internal Functions ===================

function createScatterBoxPlot( xValues, groups, ttestpValue ,figureTitle)

    % create the box plot
    boxplot( xValues, groups) 
    
    hold on
    
    % set the line width of the box plots
    set(findobj(gca,'type','line'),'linew',2)
            
    % anotate the graphs
    text(0.5, 0.9, ...
        strcat("P = ", convertPValue2SuperScript(ttestpValue ) ), ...
        'FontSize',13 ,'Unit','normalize','HorizontalAlignment','center')
    
    set(gca ,'FontSize',12, 'FontWeight','bold',...
        'LineWidth',1,'FontSize',12 ,'Box','off','XTickLabel', ...
        {'British/White','African'})
    ylabel('Levels','FontWeight','bold')
    title(figureTitle,'FontSize',13,'FontWeight','bold')
    
    % set up the colors for the plot
    colors = [0.9900 0.2250 0.0980;0.1940 0.3840 0.9560 ] ;
    
    hold on    
    % add scatter plot to the box plots
    dataSc1 = xValues(groups == false);
    dataSc2 = xValues(groups == true);
    
    % get only a smaller subset if the data points are too many
    if length(dataSc1) > 200
        dataSc1 = randsample(dataSc1,200) ;
    end
    if length(dataSc2) > 200
        dataSc2 = randsample(dataSc2,200) ;
    end
    
    % let set up the data
    x = ones(length(dataSc1)).*(1+(rand(length(dataSc1))-0.5)/5);
    x1 = ones(length(dataSc2)).*(1+(rand(length(dataSc2))-0.5)/10);
    f1 = scatter(x(:,1),dataSc1,'MarkerFaceColor', ...
        colors(1,:) , 'MarkerEdgeColor',colors(1,:) ,...
        'LineWidth',0.05);
    f1.MarkerFaceAlpha = 0.8 ;
    f1.MarkerEdgeAlpha = 0.0 ;
    hold on
    f2 = scatter(x1(:,2).*2,dataSc2,'MarkerFaceColor', ...
        colors(2,:), 'MarkerEdgeColor',colors(2,:),...
        'LineWidth',0.05);
    f2.MarkerFaceAlpha = f1.MarkerFaceAlpha;
    f2.MarkerEdgeAlpha = f1.MarkerEdgeAlpha ;
    hold off
    
end

% *********************** end of function ************************
% ****************************************************************

% ======================= another function =========================
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
    pSuperScript = sprintf('%s x 10^{%d}', firstNumbers, ...
        str2double(extractAfter(pS, 'e') )) ;

    % if the p value is large
    if p > 0.0001
        pSuperScript = sprintf('%0.4f', p) ;
    elseif p == 0
        pSuperScript = sprintf('< 1 x 10^{%d}', -300) ;
    end

end

