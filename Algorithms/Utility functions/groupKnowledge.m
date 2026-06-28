function OffspringFinal = groupKnowledge(ulOffDecs,Elite,Population,Offspring,archive1,archive2,archive3,archive4,Relate)
    Problem = PROBLEM.Current();
    OffspringFinal = [];
    for i = 1:size(ulOffDecs,1)
        ulOffDec = ulOffDecs(i,:);
        llOffDecs = repmat(Elite(unidrnd(size(Elite,2))).llDecs,Problem.Nl,1);
        dis = [];
        disE = [];
        for kE = 1:size(Elite.ulDecs,1)
            disE = [disE pdist2(ulOffDec,Elite(kE).ulDec)];
        end
        for k = 1:length(unique(Population.adds))
            dis = [dis pdist2(ulOffDec,Population((k-1)*Problem.Nl+1).ulDec)];
        end
        if ~isempty(archive2)
            % 未考虑距离
%                     if min(dis) < 0.01
            % 每个解计算距离
            for sol = 1:length(archive2)
                disE = [];
                getRelatedIndex = Relate{archive2(sol)};
                for kE = 1:size(Elite.ulDecs,1)
                    nowUseElite = Elite(kE).ulDec;
                    disE = [disE pdist2(ulOffDec(getRelatedIndex),nowUseElite(getRelatedIndex))];
                end
                [~,index] = min(disE);
                llTDec = repmat(Elite(index).llDec,Problem.Nl,1);
                llOffDecs(:,archive2(sol)) = llTDec(:,archive2(sol));
            end
%                 [~,index] = min(disE);
%                 llTDec = repmat(Elite(index).llDec,Problem.Nl,1);
%                 llOffDecs(:,archive2) = llTDec(:,archive2);
%                     end
        end
        if ~isempty(archive3)
            dis = [];
            for sol = 1:length(archive3) % 便利每个下层变量
                dis = [];
                getRelatedIndex = Relate{archive3(sol)};%获取每个下层变量相关的上层变量
                for k = 1:length(unique(Population.adds))
                    nowUsePop = Population((k-1)*Problem.Nl+1).ulDec;
                    dis = [dis pdist2(ulOffDec(getRelatedIndex),nowUsePop(getRelatedIndex))];
                end

                llTDec = Offspring((i-1)*Problem.Nl+1:i*Problem.Nl).llDecs;
                if min(dis) < 2
                    [~,index] = min(dis);
                    llTDeKno = Population((index-1)*Problem.Nl+1:index*Problem.Nl).llDecs
                    llTDec(1:Problem.Nl/2,archive3(sol)) = llTDeKno(1:Problem.Nl/2,archive3(sol));
                end
                llOffDecs(:,archive3(sol)) = llTDec(:,archive3(sol));
            end
%             llTDec = Offspring((i-1)*Problem.Nl+1:i*Problem.Nl).llDecs;
%             if min(dis) < 0.01
%                 [~,index] = min(dis);
%                 llTDeKno = Population((index-1)*Problem.Nl+1:index*Problem.Nl).llDecs
%                 llTDec(1:Problem.Nl/2) = llTDeKno(1:Problem.Nl/2);
%             end
%             llOffDecs(:,archive3) = llTDec(:,archive3);
        end
        if ~isempty(archive4)
            % 每个解计算距离
            dis = [];
            for sol = 1:length(archive4) % 便利每个下层变量
                dis = [];
                getRelatedIndex = Relate{archive4(sol)};%获取每个下层变量相关的上层变量
                for k = 1:length(unique(Population.adds))
                    nowUsePop = Population((k-1)*Problem.Nl+1).ulDec;
                    dis = [dis pdist2(ulOffDec(getRelatedIndex),nowUsePop(getRelatedIndex))];
                end

                llTDec = Offspring((i-1)*Problem.Nl+1:i*Problem.Nl).llDecs;
                if min(dis) < 2
                    [~,index] = min(dis);
                    llTDeKno = Population((index-1)*Problem.Nl+1:index*Problem.Nl).llDecs
                    llTDec(1:Problem.Nl/2,archive4(sol)) = llTDeKno(1:Problem.Nl/2,archive4(sol));
                end
                llOffDecs(:,archive4(sol)) = llTDec(:,archive4(sol));
            end

%             llTDec = Offspring((i-1)*Problem.Nl+1:i*Problem.Nl).llDecs;
%             if min(dis) < 0.01
%                 [~,index] = min(dis);
%                 llTDeKno = Population((index-1)*Problem.Nl+1:index*Problem.Nl).llDecs
%                 llTDec(1:Problem.Nl/2) = llTDeKno(1:Problem.Nl/2);
%             end
%             llOffDecs(:,archive4) = llTDec(:,archive4);
        end
        llOff = SOLUTION(ulOffDec,llOffDecs,Problem.Nu+i); 
        OffspringFinal = [OffspringFinal,llOff];
    end
end

