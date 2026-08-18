function Tasks = KnowledgeTrans(Tasks)
    Problem = PROBLEM.Current();
    % To original problem   
    transPop = cell(1,size(Tasks,2));
    for i = 2:size(Tasks,2)
        transPop{1,i} = NDSelection(Tasks(i).elitePop,'u',Problem.Nu); 
    end
    Trans = [];
    Pop1 = Tasks(1).Population;
    for i = 2:size(Tasks,2)
        transPopi = transPop{1,i};
        Popi = Tasks(i).Population;
        for j = 1:size(transPopi,2)
            ulTransDec = transPopi(j).ulDec; 
            % closest_auxiliary_elite
            elitePopi = Tasks(i).elitePop;
            tran = prod(elitePopi.ulDecs==ulTransDec,2) & prod(elitePopi.llDecs==transPopi(j).llDec,2);
            elitePopi = elitePopi(~tran);
            if size(elitePopi,2)>= Problem.Nl-1
                [~,closest] = sort(pdist2(ulTransDec,elitePopi.ulDecs));
                llTransDec = [transPopi(j).llDec; elitePopi(closest(1:Problem.Nl-1)).llDecs];
            else
                [~,closest] = sort(pdist2(ulTransDec,Popi.ulDecs));
                llTransDec = [transPopi(j).llDec; elitePopi.llDecs; Popi(closest(1:Problem.Nl-size(elitePopi,2)-1)).llDecs];
            end
            
            subTrans = SOLUTION(ulTransDec,llTransDec,Problem.Nu*(i-1)+j);
            subTrans.llEvaluate();
            subTrans = llOptimizer(subTrans);
            Tasks(1).elitePop = getElite([Tasks(1).elitePop,subTrans(subTrans.NDls==1)],size(Pop1,2));
            Trans = [Trans,subTrans];
        end
    end

    Population = [Trans,Pop1];
    nextPop = NDSelection(Population(Population.NDls==1),'u');
    nextNu = [];
    for i = 1:length(unique(nextPop.NDus))
        candiPop = nextPop(nextPop.NDus==i);
        if ~isempty(candiPop)
            if length(nextNu)< Problem.Nu && ~ismember(mode(candiPop.adds),nextNu)
                nextNu = [nextNu,mode(candiPop.adds)];
                nextPop = nextPop(nextPop.adds~=mode(candiPop.adds));
            end
            [~,rank] = sort(-candiPop.CDus);
            for k = 1:length(rank)
                if length(nextNu)< Problem.Nu  
                    if ~ismember(candiPop(rank(k)).add, nextNu)
                        nextNu = [nextNu,candiPop(rank(k)).add];
                        nextPop = nextPop(nextPop.adds~=candiPop(rank(k)).add);
                    end
                else
                    break;
                end
            end
        end
    end
    nextPop = max(Population.adds == nextNu,[],2);
    Population = Population(nextPop);
    v = 1:Problem.Nu;
    u = repelem(v,Problem.Nl);
    Population.adds(u); 
    Tasks(1).Population = Population; 
    %end original problem
end
