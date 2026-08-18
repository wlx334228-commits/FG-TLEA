function run_transfer_ratio_sensitivity()
%RUN_TRANSFER_RATIO_SENSITIVITY Run FG-TLEA Category-III/IV ratio tests.
%
% The only result data file is transfer_ratio_sensitivity_summary.csv.

    scriptPath  = mfilename('fullpath');
    projectRoot = fileparts(fileparts(scriptPath));
    cd(projectRoot);
    addpath(genpath(projectRoot));

    transferRatioValues = [0.3,0.4,0.5,0.6,0.7];
    problemNames        = {'TP1','DS1','GMP'};
    numRuns             = 10;
    interactionEpsilon  = 1e-4;

    timestamp = datestr(now,'yyyymmdd_HHMMSS');
    resultsDir = fullfile(projectRoot,'Results','transfer_ratio_sensitivity',timestamp);
    if exist(resultsDir,'dir') ~= 7
        mkdir(resultsDir);
    end

    summaryRows = {};
    writeSummary(resultsDir,summaryRows);

    cleanupParams = onCleanup(@clearExperimentParams); %#ok<NASGU>

    fprintf('FG-TLEA Category-III/IV transfer ratio sensitivity experiment\n');
    fprintf('Results directory: %s\n',resultsDir);
    fprintf('Summary: %s\n',fullfile(resultsDir,'transfer_ratio_sensitivity_summary.csv'));

    for p = 1:numel(problemNames)
        problemName = problemNames{p};
        for r = 1:numel(transferRatioValues)
            transferRatio = transferRatioValues(r);
            hvValues = nan(numRuns,1);
            runtimeValues = nan(numRuns,1);
            ulFEValues = nan(numRuns,1);
            llFEValues = nan(numRuns,1);

            for runNo = 1:numRuns
                fprintf('[%s] transferRatio=%.1f run=%d/%d\n',problemName,transferRatio,runNo,numRuns);

                global FGTLEA_PARAMS;
                FGTLEA_PARAMS = struct( ...
                    'interactionEpsilon',interactionEpsilon, ...
                    'categoryTransferRatio',transferRatio, ...
                    'skipLegacySave',true);

                rng(runNo,'twister');
                Problem = createProblem(problemName);
                Algorithm = BLEMO('outputFcn',@ignoreOutput,'save',1);

                try
                    Algorithm.Solve(Problem);
                    Population = getLastPopulation(Algorithm);
                    if ~isempty(Population)
                        hvValues(runNo) = HV(Population,Problem.optimum);
                    end
                catch ME
                    warning('FGTLEA:ExperimentFailed','%s',ME.message);
                end

                ulFEValues(runNo) = Problem.ulFE;
                llFEValues(runNo) = Problem.llFE;
                if isfield(Algorithm.metric,'runtime')
                    runtimeValues(runNo) = Algorithm.metric.runtime;
                end
            end

            validHV = ~isnan(hvValues);
            validRuntime = ~isnan(runtimeValues);
            validULFE = ~isnan(ulFEValues);
            validLLFE = ~isnan(llFEValues);

            summaryRows(end+1,:) = { ...
                problemName,transferRatio, ...
                meanOrNaN(hvValues(validHV)), ...
                stdOrNaN(hvValues(validHV)), ...
                meanOrNaN(runtimeValues(validRuntime)), ...
                meanOrNaN(ulFEValues(validULFE)), ...
                meanOrNaN(llFEValues(validLLFE))}; %#ok<AGROW>

            writeSummary(resultsDir,summaryRows);
        end
    end

    writeSummary(resultsDir,summaryRows);
    fprintf('Done. Summary: %s\n',fullfile(resultsDir,'transfer_ratio_sensitivity_summary.csv'));
end

function Problem = createProblem(problemName)
    switch problemName
        case 'TP1'
            Problem = TP1('Nu',5,'Nl',12,'Mu',2,'Ml',2,'Du',1,'Dl',2,'maxFE',500000);
        case 'DS1'
            Problem = DS1('Nu',20,'Nl',20,'Mu',2,'Ml',2,'Du',10,'Dl',10,'maxFE',1e8);
        case 'GMP'
            Problem = GMP('Nu',5,'Nl',40,'Mu',2,'Ml',2,'Du',2,'Dl',1,'maxFE',1500000);
        otherwise
            error('FGTLEA:UnknownProblem','Unknown problem: %s',problemName);
    end
end

function Population = getLastPopulation(Algorithm)
    Population = [];
    if isempty(Algorithm.result)
        return;
    end
    index = find(~cellfun(@isempty,Algorithm.result(:,2)),1,'last');
    if ~isempty(index)
        Population = Algorithm.result{index,2};
    end
end

function writeSummary(resultsDir,summaryRows)
    summaryResults = rowsToSummaryTable(summaryRows);
    writetable(summaryResults,fullfile(resultsDir,'transfer_ratio_sensitivity_summary.csv'));
end

function T = rowsToSummaryTable(rows)
    names = {'Problem','TransferRatio','HVMean','HVStd','RuntimeMean','ULFEMean','LLFEMean'};
    if isempty(rows)
        T = cell2table(cell(0,numel(names)),'VariableNames',names);
        return;
    end
    T = table(rows(:,1),cell2mat(rows(:,2)),cell2mat(rows(:,3)), ...
        cell2mat(rows(:,4)),cell2mat(rows(:,5)),cell2mat(rows(:,6)), ...
        cell2mat(rows(:,7)),'VariableNames',names);
end

function value = meanOrNaN(x)
    if isempty(x)
        value = nan;
    else
        value = mean(x);
    end
end

function value = stdOrNaN(x)
    if numel(x) < 2
        value = nan;
    else
        value = std(x);
    end
end

function ignoreOutput(~,~)
end

function clearExperimentParams()
    global FGTLEA_PARAMS;
    FGTLEA_PARAMS = [];
end
