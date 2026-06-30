function run_epsilon_sensitivity()
%RUN_EPSILON_SENSITIVITY Run FG-TLEA interaction epsilon sensitivity tests.
%
% Results are written to Results/epsilon_sensitivity/<timestamp>/, which is
% ignored by git. HV is calculated by the original project metric function.

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
    curveRows = {};
    convergenceRows = {};
    decompositionRows = {};
    writeLocationFile(resultsDir,projectRoot,epsilonValues,problemNames,numRuns);
    writeResults(resultsDir,rawRows,summaryRows,curveRows,convergenceRows,decompositionRows,epsilonValues,problemNames,numRuns);

    fprintf('FG-TLEA epsilon sensitivity experiment\n');
    fprintf('Results directory: %s\n',resultsDir);
    fprintf('Raw results: %s\n',fullfile(resultsDir,'epsilon_sensitivity_raw.csv'));
    fprintf('Summary: %s\n',fullfile(resultsDir,'epsilon_sensitivity_summary.csv'));
    fprintf('Curve: %s\n',fullfile(resultsDir,'epsilon_sensitivity_curve.csv'));
    fprintf('Convergence: %s\n',fullfile(resultsDir,'epsilon_sensitivity_convergence.csv'));
    fprintf('Decomposition: %s\n',fullfile(resultsDir,'epsilon_sensitivity_decomposition.csv'));

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

                global FGTLEA_PARAMS FGTLEA_LAST_DECOMPOSITION;
                FGTLEA_LAST_DECOMPOSITION = [];
                FGTLEA_PARAMS = struct( ...
                    'interactionEpsilon',epsilon, ...
                    'skipLegacySave',true);

                rng(runNo,'twister');
                Problem = createProblem(problemName);
                Algorithm = BLEMO('outputFcn',@recordCurve,'save',inf);
                status = 'ok';
                curveStartIndex = size(curveRows,1) + 1;
                bestHVSoFar = 0;

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

                [bestHV,fe95,fe99,auc,curvePointCount] = convergenceStats(curveRows(curveStartIndex:end,:),Problem.maxFE);
                decompositionInfo = getLastDecompositionInfo();
                convergenceRows(end+1,:) = { ...
                    problemName,epsilon,runNo,hvValue,bestHV,fe95,fe99,auc,curvePointCount,runtime,status}; %#ok<AGROW>
                decompositionRows(end+1,:) = { ...
                    problemName,epsilon,runNo,decompositionInfo.effectiveEpsilon, ...
                    decompositionInfo.du,decompositionInfo.dl, ...
                    decompositionInfo.group1Count,decompositionInfo.group2Count, ...
                    decompositionInfo.archive1Count,decompositionInfo.archive2Count, ...
                    decompositionInfo.archive3Count,decompositionInfo.archive4Count, ...
                    decompositionInfo.relateNonemptyCount,decompositionInfo.relateLinkCount, ...
                    decompositionInfo.upperMinPositive,decompositionInfo.upperMax, ...
                    decompositionInfo.lowerMinPositive,decompositionInfo.lowerMax, ...
                    decompositionInfo.status}; %#ok<AGROW>

                rawRows(end+1,:) = {problemName,epsilon,runNo,hvValue,ulFE,llFE,totalFE,runtime,status}; %#ok<AGROW>
                writeResults(resultsDir,rawRows,summaryRows,curveRows,convergenceRows,decompositionRows,epsilonValues,problemNames,numRuns);
            end

            validHV = ~isnan(hvValues);
            validFE = ~isnan(totalFEValues);
            validRuntime = ~isnan(runtimeValues);
            summaryRows(end+1,:) = { ...
                problemName,epsilon,numRuns,sum(validHV), ...
                meanOrNaN(hvValues(validHV)),stdOrNaN(hvValues(validHV)), ...
                meanOrNaN(totalFEValues(validFE)),stdOrNaN(totalFEValues(validFE)), ...
                meanOrNaN(runtimeValues(validRuntime)),stdOrNaN(runtimeValues(validRuntime))}; %#ok<AGROW>

            writeResults(resultsDir,rawRows,summaryRows,curveRows,convergenceRows,decompositionRows,epsilonValues,problemNames,numRuns);
        end
    end

    writeResults(resultsDir,rawRows,summaryRows,curveRows,convergenceRows,decompositionRows,epsilonValues,problemNames,numRuns);
    fprintf('Done. Summary: %s\n',fullfile(resultsDir,'epsilon_sensitivity_summary.csv'));

    function recordCurve(AlgorithmObj,ProblemObj)
        generation = ProblemObj.gen;
        ulFE = ProblemObj.ulFE;
        llFE = ProblemObj.llFE;
        totalFE = ulFE + llFE;

        PopulationAtGeneration = getLastPopulation(AlgorithmObj);
        curveStatus = 'ok';
        if isempty(PopulationAtGeneration)
            hvAtGeneration = nan;
            curveStatus = 'no_population';
        else
            try
                hvAtGeneration = HV(PopulationAtGeneration,ProblemObj.optimum);
                if ~isnan(hvAtGeneration)
                    bestHVSoFar = max(bestHVSoFar,hvAtGeneration);
                end
            catch ME
                hvAtGeneration = nan;
                curveStatus = ['error:',ME.identifier];
            end
        end

        runtime = nan;
        if isfield(AlgorithmObj.metric,'runtime')
            runtime = AlgorithmObj.metric.runtime;
        end

        curveRows(end+1,:) = {problemName,epsilon,runNo,generation,ulFE,llFE,totalFE,hvAtGeneration,bestHVSoFar,runtime,curveStatus}; %#ok<AGROW>
        writeResults(resultsDir,rawRows,summaryRows,curveRows,convergenceRows,decompositionRows,epsilonValues,problemNames,numRuns);
    end
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

function writeResults(resultsDir,rawRows,summaryRows,curveRows,convergenceRows,decompositionRows,epsilonValues,problemNames,numRuns)
    rawResults = rowsToRawTable(rawRows);
    summaryResults = rowsToSummaryTable(summaryRows);
    curveResults = rowsToCurveTable(curveRows);
    convergenceResults = rowsToConvergenceTable(convergenceRows);
    decompositionResults = rowsToDecompositionTable(decompositionRows);

    writetable(rawResults,fullfile(resultsDir,'epsilon_sensitivity_raw.csv'));
    writetable(summaryResults,fullfile(resultsDir,'epsilon_sensitivity_summary.csv'));
    writetable(curveResults,fullfile(resultsDir,'epsilon_sensitivity_curve.csv'));
    writetable(convergenceResults,fullfile(resultsDir,'epsilon_sensitivity_convergence.csv'));
    writetable(decompositionResults,fullfile(resultsDir,'epsilon_sensitivity_decomposition.csv'));
    save(fullfile(resultsDir,'epsilon_sensitivity_results.mat'), ...
        'rawResults','summaryResults','curveResults','convergenceResults','decompositionResults','epsilonValues','problemNames','numRuns');
end

function writeLocationFile(resultsDir,projectRoot,epsilonValues,problemNames,numRuns)
    fid = fopen(fullfile(resultsDir,'results_location.txt'),'w');
    if fid < 0
        return;
    end
    cleaner = onCleanup(@() fclose(fid));
    fprintf(fid,'Project root: %s\n',projectRoot);
    fprintf(fid,'Results directory: %s\n',resultsDir);
    fprintf(fid,'Raw CSV: %s\n',fullfile(resultsDir,'epsilon_sensitivity_raw.csv'));
    fprintf(fid,'Summary CSV: %s\n',fullfile(resultsDir,'epsilon_sensitivity_summary.csv'));
    fprintf(fid,'Curve CSV: %s\n',fullfile(resultsDir,'epsilon_sensitivity_curve.csv'));
    fprintf(fid,'Convergence CSV: %s\n',fullfile(resultsDir,'epsilon_sensitivity_convergence.csv'));
    fprintf(fid,'Decomposition CSV: %s\n',fullfile(resultsDir,'epsilon_sensitivity_decomposition.csv'));
    fprintf(fid,'MAT file: %s\n',fullfile(resultsDir,'epsilon_sensitivity_results.mat'));
    fprintf(fid,'Problems: %s\n',strjoin(problemNames,','));
    fprintf(fid,'Epsilon values: %s\n',mat2str(epsilonValues));
    fprintf(fid,'Runs per problem/epsilon: %d\n',numRuns);
    fprintf(fid,'Curve checkpoint interval: every generation\n');
    fprintf(fid,'HV metric: original Metrics/HV.m\n');
    clear cleaner;
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
    names = {'Problem','Epsilon','RunCount','ValidHVCount', ...
        'HVMean','HVStd', ...
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

function T = rowsToCurveTable(rows)
    names = {'Problem','Epsilon','Run','Generation','ULFE','LLFE','TotalFE','HV','BestHV','Runtime','Status'};
    if isempty(rows)
        T = cell2table(cell(0,numel(names)),'VariableNames',names);
        return;
    end
    T = table(rows(:,1),cell2mat(rows(:,2)),cell2mat(rows(:,3)), ...
        cell2mat(rows(:,4)),cell2mat(rows(:,5)),cell2mat(rows(:,6)), ...
        cell2mat(rows(:,7)),cell2mat(rows(:,8)),cell2mat(rows(:,9)), ...
        cell2mat(rows(:,10)),rows(:,11),'VariableNames',names);
end

function T = rowsToConvergenceTable(rows)
    names = {'Problem','Epsilon','Run','FinalHV','BestHV', ...
        'FE95Best','FE99Best','AUC','CurvePointCount','Runtime','Status'};
    if isempty(rows)
        T = cell2table(cell(0,numel(names)),'VariableNames',names);
        return;
    end
    T = table(rows(:,1),cell2mat(rows(:,2)),cell2mat(rows(:,3)), ...
        cell2mat(rows(:,4)),cell2mat(rows(:,5)),cell2mat(rows(:,6)), ...
        cell2mat(rows(:,7)),cell2mat(rows(:,8)),cell2mat(rows(:,9)), ...
        cell2mat(rows(:,10)),rows(:,11),'VariableNames',names);
end

function [bestHV,fe95,fe99,auc,curvePointCount] = convergenceStats(curveRows,maxFE)
    curvePointCount = size(curveRows,1);
    bestHV = nan;
    fe95 = nan;
    fe99 = nan;
    auc = nan;
    if isempty(curveRows)
        return;
    end

    totalFE = cell2mat(curveRows(:,7));
    hv = cell2mat(curveRows(:,9));
    valid = ~isnan(totalFE) & ~isnan(hv);
    totalFE = totalFE(valid);
    hv = hv(valid);
    if isempty(totalFE)
        return;
    end

    [totalFE,order] = sort(totalFE);
    hv = hv(order);
    bestHV = max(hv);
    if bestHV > 0
        fe95 = firstFEAtTarget(totalFE,hv,0.95*bestHV);
        fe99 = firstFEAtTarget(totalFE,hv,0.99*bestHV);
    end

    x = totalFE(:);
    y = hv(:);
    if x(1) > 0
        x = [0;x];
        y = [0;y];
    end
    if x(end) < maxFE
        x = [x;maxFE];
        y = [y;y(end)];
    end
    auc = trapz(x,y)/maxFE;
end

function fe = firstFEAtTarget(totalFE,hv,targetHV)
    index = find(hv >= targetHV,1);
    if isempty(index)
        fe = nan;
    else
        fe = totalFE(index);
    end
end

function T = rowsToDecompositionTable(rows)
    names = {'Problem','Epsilon','Run','EffectiveEpsilon','Du','Dl', ...
        'Group1Count','Group2Count','Archive1Count','Archive2Count', ...
        'Archive3Count','Archive4Count','RelateNonemptyCount','RelateLinkCount', ...
        'UpperMinPositive','UpperMax','LowerMinPositive','LowerMax','Status'};
    if isempty(rows)
        T = cell2table(cell(0,numel(names)),'VariableNames',names);
        return;
    end
    T = table(rows(:,1),cell2mat(rows(:,2)),cell2mat(rows(:,3)), ...
        cell2mat(rows(:,4)),cell2mat(rows(:,5)),cell2mat(rows(:,6)), ...
        cell2mat(rows(:,7)),cell2mat(rows(:,8)),cell2mat(rows(:,9)), ...
        cell2mat(rows(:,10)),cell2mat(rows(:,11)),cell2mat(rows(:,12)), ...
        cell2mat(rows(:,13)),cell2mat(rows(:,14)),cell2mat(rows(:,15)), ...
        cell2mat(rows(:,16)),cell2mat(rows(:,17)),cell2mat(rows(:,18)), ...
        rows(:,19),'VariableNames',names);
end

function info = getLastDecompositionInfo()
    global FGTLEA_LAST_DECOMPOSITION;
    info = defaultDecompositionInfo();
    if isempty(FGTLEA_LAST_DECOMPOSITION) || ~isstruct(FGTLEA_LAST_DECOMPOSITION)
        return;
    end
    fields = fieldnames(info);
    for i = 1:numel(fields)
        if isfield(FGTLEA_LAST_DECOMPOSITION,fields{i})
            info.(fields{i}) = FGTLEA_LAST_DECOMPOSITION.(fields{i});
        end
    end
end

function info = defaultDecompositionInfo()
    info = struct( ...
        'effectiveEpsilon',nan, ...
        'du',nan, ...
        'dl',nan, ...
        'group1Count',nan, ...
        'group2Count',nan, ...
        'archive1Count',nan, ...
        'archive2Count',nan, ...
        'archive3Count',nan, ...
        'archive4Count',nan, ...
        'relateNonemptyCount',nan, ...
        'relateLinkCount',nan, ...
        'upperMinPositive',nan, ...
        'upperMax',nan, ...
        'lowerMinPositive',nan, ...
        'lowerMax',nan, ...
        'status','missing');
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
