function elitePop = getElite(Population,N,form)
    Population = removeDupl(Population);
    if nargin<3
        form =1;
    end
    if nargin<2
        N =size(Population,2);
    end
    if form ==1 % original problem
        Population = chooseFeasiableS(Population,2);
        llElite = Population(Population.NDls == 1);
        temp = NDSelection(llElite,'u',N);
        temp = chooseFeasiableS(temp,1);
        elitePop = temp(temp.NDus == 1);
        elitePop.titles(ones(1,size(elitePop,2)));
        nTmp = temp(temp.NDus ~= 1);
        nTmp = [nTmp,Population(Population.NDls ~= 1)];
        nTmp.titles(-1*ones(1,size(nTmp,2)));
    else % bypass / many-objective
        elitePop = NDSelection(Population,'u',N);
        elitePop = elitePop(elitePop.NDus == 1);
    end
end


function Population = chooseFeasiableS(Population,form)
    i = 1;
    j = 1;
    if form == 2
        while i <= size(Population.llCons,2)
            while j <= size(Population,2) && ~isempty(Population)
                if ~isempty(Population)
                     if Population(j).llCon(i) > 0
                        Population(j) = [];
                     else
                         j = j + 1;
                     end
                end
            end
            i = i + 1;
        end
    else
        while i <= size(Population.ulCons,1)
            while j <= size(Population,2) && ~isempty(Population)
                if ~isempty(Population)
                     if Population(j).ulCon(i) > 0
                        Population(j) = [];
                     else
                         j = j+1;
                     end
                end
            end
            i = i + 1;
        end
    end
end

