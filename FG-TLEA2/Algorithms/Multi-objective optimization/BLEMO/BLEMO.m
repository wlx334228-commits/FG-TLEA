classdef BLEMO < ALGORITHM
% <multi> <real> <large/none>
% Bi-level multi-objective optimization via evolutionary multitasking
    methods
        function main(Algorithm,Problem)
            %% Generate random population
            Population = Problem.Initialization();
%             Population = llOptimizer(Population);
            Elite = getElite(Population);
            [archive1,archive2,archive3,archive4,Relate] = decomposerGroup(Elite);
            
            PPre = [];
            Pre = [];
            ulArchive = [];
            tau = 10;
            i = 0;
            ulArchive = cell(1,tau);
            
            %% Optimization
%             while Algorithm.NotTerminated(Elite,ulArchive)
%                 TmpFE = Problem.ulFE + Problem.llFE;
%                 if (TmpFE < Problem.maxFE*0.7 && mod(Problem.gen,5) ~= 0) || (TmpFE > Problem.maxFE*0.7 && mod(Problem.gen,3) ~= 0)
%                     [Population,Elite,PPre,Pre] = ulOptimizer(Population,Elite,PPre,Pre,archive1,archive2,archive3,archive4,Relate);
%                 else
%                     [Population,Elite,PPre,Pre] = ulOptimizerReal(Population,Elite,PPre,Pre);
%                 end
% 
%                 ulArchive{1,mod(i,10)+1} = Elite;
%                 i = i+1;
%             end

            while Algorithm.NotTerminated(Elite,ulArchive)
                TmpFE = Problem.ulFE + Problem.llFE;
                if Problem.gen > 10 && ((TmpFE < Problem.maxFE*0.7 && mod(Problem.gen,5) ~= 0) || (TmpFE > Problem.maxFE*0.7 && mod(Problem.gen,3) ~= 0))
                    [Population,Elite,PPre,Pre] = ulOptimizer(Population,Elite,PPre,Pre,archive1,archive2,archive3,archive4,Relate);
                else
                    [Population,Elite,PPre,Pre] = ulOptimizerReal(Population,Elite,PPre,Pre);
                end;

                ulArchive{1,mod(i,10)+1} = Elite;
                i = i+1;
            end
        end
    end
end
