function run_epsilon_sensitivity()
%RUN_EPSILON_SENSITIVITY Run FG-TLEA interaction epsilon sensitivity tests.
%
% Results are written to Results/epsilon_sensitivity/<timestamp>/, which is
% ignored by git. HV is the primary metric; FE counts and runtime are saved
% as auxiliary metrics.

    scriptPath  = mfilename('fullpath');
    projectRoot = fileparts(fileparts(scriptPath));
    cd(projectRoot);
    addpath(genpath(projectRoot));

    epsilonValues = [1e-2,1e-3,1e-4,1e-5,1e-6];
    problemNames  = {'TP1','DS1','GMP'};
    numRuns       = 15;

    timestamp = datestr(now,'yyyymmdd_HHMMSS');
    resultsDir = fullfile(projectRoot,'Results','epsilon_sensitivity',timestamp);
    if exist(resultsDir,'dir') ~= 7
        mkdir(resultsDir);
    end

    rawRows = {};
    summaryRows = {};

    fprintf('FG-TLEA epsilon sensitivity experiment\n');
    fprintf('Results directory: %s\n',resultsDir);

    for p = 1:numel(problemNames)
        problemName = problemNames{p};
        for e = 1:numel(epsilonValues)
            epsilon = epsilonValues(e);
            hvValues = nan(numRuns,1);
            ulFEValues = nan(numRuns,1);
            llFEValues = nan(numRuns,1);
            totalFEValues = nan(numRuns,1);
            runtimeValues = nan(numRuns,1);
            statuses = cell(numRuns,1);

            for runNo = 1:numRuns
                fprintf('[%s] epsilon=%g run=%d/%d\n',problemName,epsilon,runNo,numRuns);

                global FGTLEA_PARAMS;
                FGTLEA_PARAMS = struct( ...
                    'interactionEpsilon',epsilon, ...
                    'skipLegacySave',true);

                rng(runNo,'twister');
                Problem = createProblem(problemName);
                Algorithm = BLEMO('outputFcn',@silentOutput,'save',-1);
                status = 'ok';

                try
                    Algorithm.Solve(Problem);
                    Population = getLastPopulation(Algorithm);
                    if isempty(Population)
                        hvValue = nan;
                        status = 'no_population';
                    else
                        hvValue = HV(Population,Problem.optimum);
                    end
                catch ME
                    hvValue = nan;
                    status = ['error:',ME.identifier];
                    warning('FGTLEA:ExperimentFailed','%s',ME.message);
                end

                ulFE = Problem.ulFE;
                llFE = Problem.llFE;
                totalFE = ulFE + llFE;
                runtime = nan;
                if isfield(Algorithm.metric,'runtime')
                    runtime = Algorithm.metric.runtime;
                end

                hvValues(runNo) = hvValue;
                ulFEValues(runNo) = ulFE;
                llFEValues(runNo) = llFE;
                totalFEValues(runNo) = totalFE;
                runtimeValues(runNo) = runtime;
                statuses{runNo} = status;

                rawRows(end+1,:) = {problemName,epsilon,runNo,hvValue,ulFE,llFE,totalFE,runtime,status}; %#ok<AGROW>
            end

            validHV = ~isnan(hvValues);
            validFE = ~isnan(totalFEValues);
            validRuntime = ~isnan(runtimeValues);
            summaryRows(end+1,:) = { ...
                problemName,epsilon,numRuns,sum(validHV), ...
                meanOrNaN(hvValues(validHV)),stdOrNaN(hvValues(validHV)), ...
                meanOrNaN(totalFEValues(validFE)),stdOrNaN(totalFEValues(validFE)), ...
                meanOrNaN(runtimeValues(validRuntime)),stdOrNaN(runtimeValues(validRuntime))}; %#ok<AGROW>

            writeResults(resultsDir,rawRows,summaryRows,epsilonValues,problemNames,numRuns);
        end
    end

    writeResults(resultsDir,rawRows,summaryRows,epsilonValues,problemNames,numRuns);
    fprintf('Done. Summary: %s\n',fullfile(resultsDir,'epsilon_sensitivity_summary.csv'));
end

function Problem = createProblem(problemName)
    switch problemName
        case 'TP1'
            Problem = TP1('Nu',5,'Nl',12,'Mu',2,'Ml',2,'Du',1,'Dl',2,'maxFE',500000);
        case 'DS1'
            Problem = DS1('Nu',20,'Nl',20,'Mu',2,'Ml',2,'Du',10,'Dl',10,'maxFE',1500000);
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

function writeResults(resultsDir,rawRows,summaryRows,epsilonValues,problemNames,numRuns)
    rawResults = rowsToRawTable(rawRows);
    summaryResults = rowsToSummaryTable(summaryRows);

    writetable(rawResults,fullfile(resultsDir,'epsilon_sensitivity_raw.csv'));
    writetable(summaryResults,fullfile(resultsDir,'epsilon_sensitivity_summary.csv'));
    save(fullfile(resultsDir,'epsilon_sensitivity_results.mat'), ...
        'rawResults','summaryResults','epsilonValues','problemNames','numRuns');
end

function T = rowsToRawTable(rows)
    names = {'Problem','Epsilon','Run','HV','ULFE','LLFE','TotalFE','Runtime','Status'};
    if isempty(rows)
        T = cell2table(cell(0,numel(names)),'VariableNames',names);
        return;
    end
    T = table(rows(:,1),cell2mat(rows(:,2)),cell2mat(rows(:,3)), ...
        cell2mat(rows(:,4)),cell2mat(rows(:,5)),cell2mat(rows(:,6)), ...
        cell2mat(rows(:,7)),cell2mat(rows(:,8)),rows(:,9), ...
        'VariableNames',names);
end

function T = rowsToSummaryTable(rows)
    names = {'Problem','Epsilon','RunCount','ValidHVCount','HVMean','HVStd', ...
        'TotalFEMean','TotalFEStd','RuntimeMean','RuntimeStd'};
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

function silentOutput(~,~)
end
