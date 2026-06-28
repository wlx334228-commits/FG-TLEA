function [Population,Elite,PPre,Pre] = ulOptimizer(Population,Elite,PPre,Pre,archive1,archive2,archive3,archive4,Relate)
    Problem = PROBLEM.Current();
    Offspring = [];
    ulOffDecs = [];
    for i = 1:Problem.Nu
        subP = Population(Population.adds==i);
        DecP1 = subP(1).ulDec;
        eliteP = Elite(unidrnd(size(Elite,2)));
        DecP2 = eliteP.ulDec;
        restP = Population(Population.adds~=i);  
        DecP3 = restP(unidrnd(size(restP,2))).ulDec;            
        ulOffDec = OperatorDE(DecP1,DecP2,DecP3,'u');
        
        %% 使用DE差分进化生成下层变量
        Decpl1 = subP.llDecs;
        eli = Elite(unidrnd(size(Elite,2))).llDecs;
        Decpl2 = repmat(eli,size(Decpl1,1),1);
        while 1
            randP = unidrnd(Problem.Nu);
            if randP ~= i
                Decpl3 = Population((randP-1)*Problem.Nl+1:randP*Problem.Nl).llDecs;
                break;
            end
        end
        llOffDecs = OperatorDE(Decpl1,Decpl2,Decpl3,'l');
        
        
        %% 使用PCX + PM生成下层变量
%         llOffDecs = selec_recomb_muta(Population,Problem.Nl,Problem.llMin,Problem.llMax);
% 
        llOff = SOLUTION(ulOffDec,llOffDecs,Problem.Nu+i);

        ulOffDecs = [ulOffDecs;ulOffDec];
        Offspring = [Offspring,llOff];
%         if Problem.gen >= 3 && (mod(Problem.gen,3) == 0  || isequal(sort(unique(PPre.ulDecs,"rows")),sort(unique(Pre.ulDecs,"rows"))))
% %             PPre = Pre;
% %             Pre = Population;
%             llOffDecs = selec_recomb_muta(Population,Problem.Nl,Problem.llMin,Problem.llMax);
%             llOff = SOLUTION(ulOffDec,llOffDecs,Problem.Nu+i); 
% %             llOff.llEvaluate();
% %             tic;
% %             llOff = llOptimizern(llOff);
% %             toc;
% %             Elite = getElite([Elite,llOff],Problem.Nu*Problem.Nl);
%             Offspring = [Offspring,llOff];
% %             [Offspring,Elite] = multi_task(Offspring,Elite,Population);
%         else
% %         llOffDecs = selec_recomb_muta(Population,Problem.Nl,Problem.llMin,Problem.llMax);
% %         llOff = SOLUTION(ulOffDec,llOffDecs,Problem.Nu+i); 
% %         llOff.llEvaluate();
% %         tic;
% %         llOff = llOptimizer(llOff);
% %         toc;
% %         Elite = getElite([Elite,llOff],Problem.Nu*Problem.Nl);
% %         Offspring = [Offspring,llOff];
%             ulOffDecs = [ulOffDecs;ulOffDec];
%         end
    end

    %% surrogate model
    %train the surrogate model
    if Problem.gen < 3 || (mod(Problem.gen,3) ~= 0 && ~isequal(sort(unique(PPre.ulDecs,"rows")),sort(unique(Pre.ulDecs,"rows"))))
        if Problem.gen == 1
            PPre = Population;
            Pre = Population;
        elseif Problem.gen == 2
            Pre = Population;
        else
            PPre = Pre;
            Pre = Population;
        end
        pos = Elite;
        trainlabel = ones(size(unique(Elite.ulDecs,"rows"),1),1);
        traindata = unique(Elite.ulDecs,"rows");
        neg = setdiff(unique(Population.ulDecs,"rows"),unique(Elite.ulDecs,"rows"),"rows");
        trainlabel = [trainlabel;-1 * ones(size(neg,1),1)];
        traindata = [traindata;neg];
        model = train(trainlabel,sparse(traindata));
    
        reslabel = ulOffDecs * model.w';
        if sum(reslabel > 0,1) >= ceil(Problem.Nu*0.4)
            A = find(reslabel > 0);
            ot = [];
            for i = 1:size(A,1)
                ot = [ot Offspring((A(i)-1)*Problem.Nl+1:A(i)*Problem.Nl)];
            end
            Offspring = ot;
            ulOffDecs = ulOffDecs(reslabel > 0,:);
        else
            [~,A] = sort(reslabel,"descend");
            ot = [];
            for i = 1:ceil(Problem.Nu*0.4)
                ot = [ot Offspring((A(i)-1)*Problem.Nl+1:A(i)*Problem.Nl)];
            end
            Offspring = ot;
            ulOffDecs = ulOffDecs(A(1:ceil(Problem.Nu*0.4)),:);
        end
        
        % 根据分组学习
        Offspring = groupKnowledge(ulOffDecs,Elite,Population,Offspring,archive1,archive2,archive3,archive4,Relate);
        [Offspring,Elite,calHvS] = multi_task(Offspring,Elite,Population,archive1,archive2);
        calHV = HV(calHvS,Problem.optimum);
    else
        PPre = Pre;
        Pre = Population;
        Offspring = groupKnowledge(ulOffDecs,Elite,Population,Offspring,archive1,archive2,archive3,archive4,Relate);
        [Offspring,Elite,calHvS] = multi_task(Offspring,Elite,Population,archive1,archive2);
        calHV = HV(calHvS,Problem.optimum);
    end

    %%
    namestr = ['h' num2str(Problem.gen)  '=calHV'];
    eval(namestr);
    if ~FGTLEAGetParameter('skipLegacySave',false)
        if isa(Problem,'TP1')
            save('TP1HV.mat', ['h' num2str(Problem.gen)], '-append');
        elseif isa(Problem,'TP2')
            save('TP2HV.mat', ['h' num2str(Problem.gen)], '-append');
        elseif isa(Problem,'DS1')
            save('DS1HV.mat', ['h' num2str(Problem.gen)], '-append');
        elseif isa(Problem,'DS2')
            save('DS2HV.mat', ['h' num2str(Problem.gen)], '-append');
        elseif isa(Problem,'DS3')
            save('DS3HV.mat', ['h' num2str(Problem.gen)], '-append');
        elseif isa(Problem,'DS4')
            save('DS4HV.mat', ['h' num2str(Problem.gen)], '-append');
        elseif isa(Problem,'DS5')
            save('DS5HV.mat', ['h' num2str(Problem.gen)], '-append');
        end
    end

    FE = Problem.ulFE + Problem.llFE;
    namestr = ['f' num2str(Problem.gen) '=FE'];
    eval(namestr);
    if ~FGTLEAGetParameter('skipLegacySave',false)
        if isa(Problem,'TP1')
            save('TP1FE.mat', ['f' num2str(Problem.gen)], '-append');
        elseif isa(Problem,'TP2')
            save('TP2FE.mat', ['f' num2str(Problem.gen)], '-append');
        elseif isa(Problem,'DS1')
            save('DS1FE.mat', ['f' num2str(Problem.gen)], '-append');
        elseif isa(Problem,'DS2')
            save('DS2FE.mat', ['f' num2str(Problem.gen)], '-append');
        elseif isa(Problem,'DS3')
            save('DS3FE.mat', ['f' num2str(Problem.gen)], '-append');
        elseif isa(Problem,'DS4')
            save('DS4FE.mat', ['f' num2str(Problem.gen)], '-append');
        elseif isa(Problem,'DS5')
            save('DS5FE.mat', ['f' num2str(Problem.gen)], '-append');
        end
    end

    %%
    Population = [Population,Offspring];
    nextPop = NDSelection(Population(Population.NDls==1),'u');
    eTmp = nextPop(nextPop.NDus==1);
    eTmp.titles(ones(1,size(eTmp,2)));
    nTmp = nextPop(nextPop.NDus ~= 1);
    nTmp.titles(-1*ones(1,size(nTmp,2)));
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
end
function Decs = selec_recomb_muta(Population,N,Lower,Upper)
    MatingPool = TournamentSelection(2,N,Population.NDls,-Population.CDls);
    Parent = Population(MatingPool).llDecs;
    [proC,proM,disM] = deal(1,1,20);
    [N,D] = size(Parent);

    %% PCX
    W         = mean(Parent,1);
    p1        = [2:N,1];
    p2        = [3:N,1,2];
    Offspring = Parent + 0.1*(Parent-W) + mean(abs(Parent-W),2).*(Parent(p2,:)-Parent(p1,:))/2;
    Site      = repmat(rand(N,1)>proC,1,D);
    Offspring(Site) = Parent(Site);
    
    %% Polynomial mutation
    Lower = repmat(Lower,N,1);
    Upper = repmat(Upper,N,1);
    Site  = rand(N,D) < proM/D;
    mu    = rand(N,D);
    temp  = Site & mu<=0.5;
    Offspring       = min(max(Offspring,Lower),Upper);
    Offspring(temp) = Offspring(temp)+(Upper(temp)-Lower(temp)).*((2.*mu(temp)+(1-2.*mu(temp)).*...
                      (1-(Offspring(temp)-Lower(temp))./(Upper(temp)-Lower(temp))).^(disM+1)).^(1/(disM+1))-1);
    temp = Site & mu>0.5; 
    Offspring(temp) = Offspring(temp)+(Upper(temp)-Lower(temp)).*(1-(2.*(1-mu(temp))+2.*(mu(temp)-0.5).*...
                      (1-(Upper(temp)-Offspring(temp))./(Upper(temp)-Lower(temp))).^(disM+1)).^(1/(disM+1)));
    Decs = Offspring;
end

function [Offspring,Elite,calHv] = multi_task(Offspring,Elite,Population,archive1,archive2)
    if isempty(Offspring)
        calHv = 0;
        return;
    end
    Problem = PROBLEM.Current();
    rmp=1; % Random mating probability
    gen = 250; % Maximum Number of generations
    muc = 10; % Distribution Index of SBX crossover operator
    mum = 10; % Distribution Index of Polynomial Mutation operator
    prob_vswap = 0;
    
    %% MO-MFEA
    % divide several groups -- one group can over 2 task
    % groups = (1:Problem.Nu);
    % for i = 1:1 
    %     if ~ismember(i,groups)
    %         continue
    %     end
    %     distance = -ones(1,Problem.Nu);
    %     for j = 1:Problem.Nu
    %         if i ~= j
    %             distance(j) = sum((Offspring(1,12*(i-1)+1).ulDec - Offspring(1,12*(j-1)+1).ulDec).*(Offspring(1,12*(i-1)+1).ulDec - Offspring(1,12*(j-1)+1).ulDec),2);  
    %         end
    %     end
    % 
    %     group{i} = find(distance < 0.1);
    %     groups = groups(~group{i});
    % end
    
    % divide several groups -- one group only 2 task
    len = size(unique(Offspring.ulDecs,"rows"),1);
    groups = (1:len);
    group = {};
    tmp = 1;
    while ~isempty(groups)
        if length(groups) == 1
            group{tmp} = groups(1);
            break
        end
        distance = [];
        for i = 2:length(groups)
            distance(i-1) = sum((Offspring(1,12*(groups(1)-1)+1).ulDec - Offspring(1,12*(groups(i)-1)+1).ulDec).*(Offspring(1,12*(groups(1)-1)+1).ulDec - Offspring(1,12*(groups(i)-1)+1).ulDec),2);  
        end
        if min(distance) <= 0.005
            [~,minIndex] = min(distance);
            group{tmp} = [groups(1) groups(minIndex+1)];
            groups(minIndex+1) = [];
        else
            group{tmp} = [groups(1)];
        end
%         [~,minIndex] = min(distance);
%         group{tmp} = [groups(1) groups(minIndex+1)];
%         groups(minIndex+1) = [];
        groups(1) = [];
        tmp = tmp + 1;
    end

    % 多任务优化
    for i = 1:length(group)
        if size(group{i},2) == 1
            beginIndex = (group{i}(1)-1)*Problem.Nl;
            llOff = Offspring((beginIndex+1):(beginIndex+Problem.Nl));
            llOff.llEvaluate();
            llOff = llOptimizerArchive(llOff,archive1,archive2);
            llOff.titles(-1*ones(size(llOff,2),1));
            Offspring((beginIndex+1):(beginIndex+Problem.Nl)) = llOff;
            calHv = getNewElite([Elite,llOff],Problem.Nu*Problem.Nl);
            Elite = getElite([Elite,llOff],Problem.Nu*Problem.Nl);
            continue;
        end
        G = 0;
        llArchive1 = [];
        llArchive2 = [];
        tau = 10;

        beginIndex1 = (group{i}(1)-1)*Problem.Nl;
        population_T1 = Offspring((beginIndex1+1):(beginIndex1+Problem.Nl));
        beginIndex2 = (group{i}(2)-1)*Problem.Nl;
        population_T2 = Offspring((beginIndex2+1):(beginIndex2+Problem.Nl));
        population_T1.llEvaluate();
        population_T1 = NDSelection(population_T1,'l',Problem.Nl);
        population_T2.llEvaluate();
        population_T2 = NDSelection(population_T2,'l',Problem.Nl);

        population = [population_T1 population_T2];
        pop = length(population);
        dim = max(length(population_T1(1).llDec),length(population_T2(1).llDec));
        
        while terminationCheck(llArchive1,G) && terminationCheck(llArchive2,G)
            llArchive1 = cell(1,tau);
            llArchive2 = cell(1,tau);
            
            for t = 1:tau
                MatingPoll = TournamentSelection(2,size(population,2),population.rankl);
                parent = population(MatingPoll);
%                 for i = 1:pop % Performing binary tournament selection to create parent pool
%                     p1=1+round(rand(1)*(pop-1));
%                     p2=1+round(rand(1)*(pop-1));
%                     if population(p1).rankl < population(p2).rankl
%                         parent(i) = SOLUTION(population(p1).ulDec,population(p1).llDec,population(p1).add);
%                     elseif population(p1).rankl == population(p2).rankl
%                         if rand(1) <= 0.5
%                             parent(i) = SOLUTION(population(p1).ulDec,population(p1).llDec,population(p1).add);
%                         else
%                             parent(i) = SOLUTION(population(p2).ulDec,population(p2).llDec,population(p2).add);
%                         end
%                     else
%                        parent(i) = SOLUTION(population(p2).ulDec,population(p2).llDec,population(p2).add);
%                     end
%                 end
    
                count=1;
                for i=1:2:pop-1 % Create offspring population via mutation and crossover
                    child(count)=SOLUTION();
                    child(count+1)=SOLUTION();
                    p1=i;
                    p2=i+1;
                    minAdd = min(parent(p1).add,parent(p2).add);
                    maxAdd = max(parent(p1).add,parent(p2).add);
                    if parent(p1).add==parent(p2).add
                        % % DE
                        parentOff = Offspring(Offspring.adds == parent(p1).add);
                        elitePop = parentOff(parentOff.NDls==1);
                        Decp2_1 = elitePop((unidrnd(size(elitePop,2)))).llDec;
                        while 1
                            randP = unidrnd(Problem.Nl);
                            if parentOff(randP) ~= parent(p1)
                                Decp3_1 = parentOff(randP).llDec;
                                break;
                            end
                        end
                        child(count).llDec = OperatorDE(parent(p1).llDec,Decp2_1,Decp3_1,'l');
                        child(count).llDec([archive1,archive2]) = parent(p1).llDec([archive1,archive2]);
                        child(count).add = parent(p1).add;
                        child(count).ulDec = parent(p1).ulDec;
                        
                        Decp2_2 = elitePop((unidrnd(size(elitePop,2)))).llDec;
                        while 1
                            randP = unidrnd(Problem.Nl);
                            if parentOff(randP) ~= parent(p1)
                                Decp3_2 = parentOff(randP).llDec;
                                break;
                            end
                        end
                        child(count+1).llDec = OperatorDE(parent(p2).llDec,Decp2_2,Decp3_2,'l');
                        child(count+1).llDec([archive1,archive2]) = parent(p2).llDec([archive1,archive2]);
                        child(count+1).add = parent(p2).add;
                        child(count+1).ulDec = parent(p2).ulDec;


                        % SBX PM
%                         [child(count).llDec,child(count+1).llDec]=crossover(parent(p1).llDec,parent(p2).llDec,muc,dim,prob_vswap);
%                         child(count).llDec = mutate(child(count).llDec,mum,dim,1/dim);
%                         child(count+1).llDec=mutate(child(count+1).llDec,mum,dim,1/dim);
%                         child(count).add=parent(p1).add;
%                         child(count).ulDec=parent(p1).ulDec;
%                         child(count+1).add=parent(p2).add;
%                         child(count+1).ulDec=parent(p2).ulDec;
                    else
                        if rand(1)<rmp
                            %DE
                            parentOff1 = Offspring(Offspring.adds == parent(p1).add);
                            parentOff2 = Offspring(Offspring.adds == parent(p2).add);
                            elitePop1 = parentOff1(parentOff1.NDls==1);
                            elitePop2 = parentOff2(parentOff2.NDls==1);
                            Decpt1_2 = elitePop1((unidrnd(size(elitePop1,2)))).llDec;
                            Decpt2_2 = elitePop2((unidrnd(size(elitePop2,2)))).llDec;
                            
                            while 1
                                randP = unidrnd(Problem.Nl);
                                if parentOff1(randP) ~= parent(p1)
                                    Decpt1_3 = parentOff(randP).llDec;
                                    break;
                                end
                            end
                            while 1
                                randP = unidrnd(Problem.Nl);
                                if parentOff2(randP) ~= parent(p2)
                                    Decpt2_3 = parentOff(randP).llDec;
                                    break;
                                end
                            end
%                             randP1 = unidrnd(Problem.Nl);
%                             Decpt1_3 = parentOff1(randP1).llDec;
%                             
%                             randP2 = unidrnd(Problem.Nl);
%                             Decpt2_3 = parentOff2(randP2).llDec;
                            
                            child(count).llDec = OperatorDE(parent(p1).llDec,Decpt2_2,Decpt1_3,'l');
                            child(count).llDec([archive1,archive2]) = parent(p1).llDec([archive1,archive2]);

                            child(count).add = parent(p1).add;
                            child(count).ulDec = parent(p1).ulDec;
                            child(count+1).llDec = OperatorDE(parent(p2).llDec,Decpt1_2,Decpt2_3,'l');
                            child(count+1).llDec([archive1,archive2]) = parent(p2).llDec([archive1,archive2]);
                            child(count+1).add = parent(p2).add;
                            child(count+1).ulDec = parent(p2).ulDec;
   
                            
                            %SBX PM
%                             [child(count).llDec,child(count+1).llDec]= crossover(parent(p1).llDec,parent(p2).llDec,muc,dim,prob_vswap); 
%                             child(count).llDec = mutate(child(count).llDec,mum,dim,1/dim);
%                             child(count+1).llDec=mutate(child(count+1).llDec,mum,dim,1/dim);
%                             child(count).add=round(rand(1))*(maxAdd-minAdd)+minAdd;
%                             anst = parent(parent.adds==child(count).add);
%                             child(count).ulDec=anst(1).ulDec;
%                             child(count+1).add=round(rand(1))*(maxAdd-minAdd)+minAdd;
%                             anst = parent(parent.adds==child(count+1).add);
%                             child(count+1).ulDec=anst(1).ulDec;
                        else
                            child(count).llDec = mutate(parent(p1).llDec,mum,dim,1);
                            child(count).llDec([archive1,archive2]) = parent(p1).llDec([archive1,archive2]);
                            child(count+1).llDec=mutate(parent(p2).llDec,mum,dim,1);
                            child(count+1).llDec([archive1,archive2]) = parent(p2).llDec([archive1,archive2]);
                            child(count).add=parent(p1).add;
                            child(count).ulDec=parent(p1).ulDec;
                            child(count+1).add=parent(p2).add;
                            child(count+1).ulDec=parent(p2).ulDec;
                        end
                    end
                    count=count+2;
                end
                child.llEvaluate();
                %population=reset(population,pop);
                intpopulation(1:pop)=population;
                intpopulation(pop+1:2*pop)=child;
                intpopulation_T1=intpopulation([intpopulation.add]==population_T1(1).add);
                intpopulation_T2=intpopulation([intpopulation.add]==population_T2(1).add);
                T1_pop=length(intpopulation_T1);
                T2_pop=length(intpopulation_T2);
                intpopulation_T1 = NDSelection(intpopulation_T1,'l',Problem.Nl);
                intpopulation_T2 = NDSelection(intpopulation_T2,'l',Problem.Nl);
                population = [intpopulation_T1 intpopulation_T2];
                Offspring((beginIndex1+1):(beginIndex1+Problem.Nl)) = intpopulation_T1;
                elitePop1 = intpopulation_T1(intpopulation_T1.NDls==1);
                llArchive1{1,t} = elitePop1;
                
                Offspring((beginIndex2+1):(beginIndex2+Problem.Nl)) = intpopulation_T2;
                elitePop2 = intpopulation_T2(intpopulation_T2.NDls==1);
                llArchive2{1,t} = elitePop2;
                
                G=G+1;
            end   
        end 
        t1OffDecs = selec_recomb_muta(Population,Problem.Nl,Problem.llMin,Problem.llMax);
        t1Off = SOLUTION(unique(intpopulation_T1.ulDecs,"rows"),t1OffDecs,unique(intpopulation_T1.adds));
        t1 = [t1Off intpopulation_T1];
        t1.llEvaluate();
        idx1 = randperm(numel(t1));
        t2OffDecs = selec_recomb_muta(Population,Problem.Nl,Problem.llMin,Problem.llMax);
        t2Off = SOLUTION(unique(intpopulation_T2.ulDecs,"rows"),t2OffDecs,unique(intpopulation_T2.adds));
        t2 = [t2Off intpopulation_T2];
        t2.llEvaluate();
        idx2 = randperm(numel(t2));
        intpopulation_T1 = llOptimizerArchive(t1(idx1(1:Problem.Nl)),archive1,archive2);
        intpopulation_T2 = llOptimizerArchive(t2(idx2(1:Problem.Nl)),archive1,archive2);
%         Offspring((beginIndex1+1):(beginIndex1+Problem.Nl)) = intpopulation_T1;
%         Offspring((beginIndex2+1):(beginIndex2+Problem.Nl)) = intpopulation_T2;
        intpopulation_T1.titles(-1*ones(size(intpopulation_T1,2),1));
        intpopulation_T2.titles(-1*ones(size(intpopulation_T2,2),1));
        calHv = getNewElite([Elite,intpopulation_T1,intpopulation_T2],Problem.Nu*Problem.Nl);
        Elite = getElite([Elite,intpopulation_T1,intpopulation_T2],Problem.Nu*Problem.Nl);
        Offspring((beginIndex1+1):(beginIndex1+Problem.Nl)) = intpopulation_T1;
        Offspring((beginIndex2+1):(beginIndex2+Problem.Nl)) = intpopulation_T2;
    end 
end
