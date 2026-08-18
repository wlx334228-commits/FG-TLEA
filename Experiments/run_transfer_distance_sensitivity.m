function run_transfer_distance_sensitivity()
%RUN_TRANSFER_DISTANCE_SENSITIVITY Run FG-TLEA transfer-distance tests.
%
% Results include one raw CSV for each run and one summary CSV.

    scriptPath  = mfilename('fullpath');
    projectRoot = fileparts(fileparts(scriptPath));
    cd(projectRoot);
    addpath(genpath(projectRoot));

    distanceThresholdValues = [0.005,0.01,0.015,0.2];
    problemNames            = {'DS1'};
    numRuns                 = 4;
    interactionEpsilon      = 1e-4;
    transferRatio           = 0.5;

    timestamp = datestr(now,'yyyymmdd_HHMMSS');
    resultsDir = fullfile(projectRoot,'Results','transfer_distance_sensitivity',timestamp);
    if exist(resultsDir,'dir') ~= 7
        mkdir(resultsDir);
    end

    rawRows = {};
    summaryRows = {};
    writeResults(resultsDir,rawRows,summaryRows);

    cleanupParams = onCleanup(@clearExperimentParams); %#ok<NASGU>

    fprintf('FG-TLEA Category-III/IV transfer distance sensitivity experiment\n');
    fprintf('Results directory: %s\n',resultsDir);
    fprintf('Raw results: %s\n',fullfile(resultsDir,'transfer_distance_sensitivity_raw.csv'));
    fprintf('Summary: %s\n',fullfile(resultsDir,'transfer_distance_sensitivity_summary.csv'));

    for p = 1:numel(problemNames)
        problemName = problemNames{p};
        for d = 1:numel(distanceThresholdValues)
            distanceThreshold = distanceThresholdValues(d);
            hvValues = nan(numRuns,1);
            runtimeValues = nan(numRuns,1);
            ulFEValues = nan(numRuns,1);
            llFEValues = nan(numRuns,1);
            attemptCounts = nan(numRuns,1);
            acceptedCounts = nan(numRuns,1);
            rejectedCounts = nan(numRuns,1);

            for runNo = 1:numRuns
                fprintf('[%s] distanceThreshold=%g run=%d/%d\n',problemName,distanceThreshold,runNo,numRuns);

                global FGTLEA_PARAMS FGTLEA_TRANSFER_DIAGNOSTICS;
                FGTLEA_TRANSFER_DIAGNOSTICS = defaultTransferDiagnostics();
                FGTLEA_PARAMS = struct( ...
                    'interactionEpsilon',interactionEpsilon, ...
                    'categoryTransferRatio',transferRatio, ...
                    'categoryTransferDistanceThreshold',distanceThreshold, ...
                    'trackTransferDiagnostics',true, ...
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
                diagnostics = getTransferDiagnostics();
                attemptCounts(runNo) = diagnostics.attemptCount;
                acceptedCounts(runNo) = diagnostics.acceptedCount;
                rejectedCounts(runNo) = diagnostics.rejectedCount;

                rawRows(end+1,:) = { ...
                    problemName,distanceThreshold,runNo,hvValues(runNo), ...
                    runtimeValues(runNo),ulFEValues(runNo),llFEValues(runNo), ...
                    diagnostics.attemptCount,diagnostics.acceptedCount,diagnostics.rejectedCount}; %#ok<AGROW>
                writeResults(resultsDir,rawRows,summaryRows);
            end

            validHV = ~isnan(hvValues);
            validRuntime = ~isnan(runtimeValues);
            validULFE = ~isnan(ulFEValues);
            validLLFE = ~isnan(llFEValues);
            validAttempts = ~isnan(attemptCounts);
            validAccepted = ~isnan(acceptedCounts);
            validRejected = ~isnan(rejectedCounts);

            summaryRows(end+1,:) = { ...
                problemName,distanceThreshold, ...
                meanOrNaN(hvValues(validHV)), ...
                stdOrNaN(hvValues(validHV)), ...
                meanOrNaN(runtimeValues(validRuntime)), ...
                meanOrNaN(ulFEValues(validULFE)), ...
                meanOrNaN(llFEValues(validLLFE)), ...
                meanOrNaN(attemptCounts(validAttempts)), ...
                meanOrNaN(acceptedCounts(validAccepted)), ...
                meanOrNaN(rejectedCounts(validRejected))}; %#ok<AGROW>

            writeResults(resultsDir,rawRows,summaryRows);
        end
    end

    writeResults(resultsDir,rawRows,summaryRows);
    fprintf('Done. Summary: %s\n',fullfile(resultsDir,'transfer_distance_sensitivity_summary.csv'));
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

function writeResults(resultsDir,rawRows,summaryRows)
    rawResults = rowsToRawTable(rawRows);
    summaryResults = rowsToSummaryTable(summaryRows);
    writetable(rawResults,fullfile(resultsDir,'transfer_distance_sensitivity_raw.csv'));
    writetable(summaryResults,fullfile(resultsDir,'transfer_distance_sensitivity_summary.csv'));
end

function T = rowsToRawTable(rows)
    names = {'Problem','DistanceThreshold','Run','HV','Runtime','ULFE','LLFE', ...
        'TransferAttemptCount','TransferAcceptedCount','TransferRejectedCount'};
    if isempty(rows)
        T = cell2table(cell(0,numel(names)),'VariableNames',names);
        return;
    end
    T = table(rows(:,1),cell2mat(rows(:,2)),cell2mat(rows(:,3)), ...
        cell2mat(rows(:,4)),cell2mat(rows(:,5)),cell2mat(rows(:,6)), ...
        cell2mat(rows(:,7)),cell2mat(rows(:,8)),cell2mat(rows(:,9)), ...
        cell2mat(rows(:,10)),'VariableNames',names);
end

function T = rowsToSummaryTable(rows)
    names = {'Problem','DistanceThreshold','HVMean','HVStd','RuntimeMean','ULFEMean','LLFEMean', ...
        'TransferAttemptCountMean','TransferAcceptedCountMean','TransferRejectedCountMean'};
    if isempty(rows)
        T = cell2table(cell(0,numel(names)),'VariableNames',names);
        return;
    end
    T = table(rows(:,1),cell2mat(rows(:,2)),cell2mat(rows(:,3)), ...
        cell2mat(rows(:,4)),cell2mat(rows(:,5)),cell2mat(rows(:,6)), ...
        cell2mat(rows(:,7)),cell2mat(rows(:,8)),cell2mat(rows(:,9)), ...
        cell2mat(rows(:,10)),'VariableNames',names);
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

function diagnostics = getTransferDiagnostics()
    global FGTLEA_TRANSFER_DIAGNOSTICS;
    diagnostics = defaultTransferDiagnostics();
    if isempty(FGTLEA_TRANSFER_DIAGNOSTICS) || ~isstruct(FGTLEA_TRANSFER_DIAGNOSTICS)
        return;
    end

    fields = fieldnames(diagnostics);
    for i = 1:numel(fields)
        if isfield(FGTLEA_TRANSFER_DIAGNOSTICS,fields{i})
            diagnostics.(fields{i}) = FGTLEA_TRANSFER_DIAGNOSTICS.(fields{i});
        end
    end
end

function diagnostics = defaultTransferDiagnostics()
    diagnostics = struct('attemptCount',0,'acceptedCount',0,'rejectedCount',0);
end

function clearExperimentParams()
    global FGTLEA_PARAMS FGTLEA_TRANSFER_DIAGNOSTICS;
    FGTLEA_PARAMS = [];
    FGTLEA_TRANSFER_DIAGNOSTICS = [];
end
