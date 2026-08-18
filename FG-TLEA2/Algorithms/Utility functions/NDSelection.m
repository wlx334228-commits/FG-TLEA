function Population = NDSelection(Population,flag,N)
% The environmental selection of NSGA-II

%------------------------------- Copyright --------------------------------
% Copyright (c) 2021 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------
    if nargin<3
        N = size(Population,2);
    else
        N = min(size(Population,2),N);
    end
    if flag == 'u'
        for i = 1:size(Population,2)
            if isempty(Population(i).ulObj)
                Population(i).ulEvaluate();
            end
        end
        Objs = Population.ulObjs;
        Con = Population.ulCons;
    elseif flag == 'l'
        for i = 1:size(Population,2)
            if isempty(Population(i).llObj)
                Population(i).llEvaluate();
            end
        end
        Objs = Population.llObjs;
        Con = Population.llCons;
    end
    %% Non-dominated sorting
    [FrontNo,MaxFNo] = NDSort(Objs,Con,N);
    Next = FrontNo < MaxFNo;
    
    %% Calculate the crowding distance of each solution
    CrowdDis = CrowdingDistance(Objs,FrontNo);

    %% additional code
    rank = 1;
    for i = 1:MaxFNo
        index = find(FrontNo == i);
        [~,reflect] = sort(CrowdDis(index),2,"descend");
        index = index(reflect);
        for j = 1 : length(index)
            Population(1,index(j)).rankl = rank;
            rank = rank + 1;
        end
    end
    
    %% Select the solutions in the last front based on their crowding distances
    Last     = find(FrontNo==MaxFNo);
    [~,Rank] = sort(CrowdDis(Last),'descend');
    Next(Last(Rank(1:N-sum(Next)))) = true;
    
    %% Population for next generation
    Population = Population(Next);
    FrontNo    = FrontNo(Next);
    CrowdDis   = CrowdDis(Next);

    if flag == 'u'
        Population.NDus(FrontNo');
        Population.CDus(CrowdDis');
    elseif flag == 'l'
        Population.NDls(FrontNo');
        Population.CDls(CrowdDis');
    end
end