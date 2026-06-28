function [archive1,archive2,archive3,archive4,Relate] = decomposer(basicSolution)
    Problem = PROBLEM.Current;
    basicsolution = basicSolution(unidrnd(size(basicSolution,2)));
    ulDec = basicsolution.ulDec;
    llDec = basicsolution.llDec;
    epsilon = FGTLEAGetParameter('interactionEpsilon',1e-4);
    Dec = [basicsolution.ulDec basicsolution.llDec];

    % 将上层每一维决策变量进行改变
    upperVariedSolution = [];
    for i = 1:Problem.Du
        s = ulDec;
        s(i) = s(i)*1.5;
        us = SOLUTION(s,llDec);
        us.llEvaluate();
        us.ulEvaluate();
        upperVariedSolution = [upperVariedSolution us];
    end

    lowerVariedSolution = [];
    for i = 1:Problem.Dl
        s = llDec;
        s(i) = s(i)*1.5;
        ls = SOLUTION(ulDec,s);
        ls.llEvaluate();
        ls.ulEvaluate();
        lowerVariedSolution = [lowerVariedSolution ls];
    end
    
    bothVariedSolution = [];
    for i = 1:Problem.Du
        for j = 1:Problem.Dl
            s = Dec;
            s(i) = s(i)*1.5;
            s(Problem.Du+j) = s(Problem.Du+j)*1.5;
            bs = SOLUTION(s(1:Problem.Du),s(Problem.Du+1:end));
            bs.llEvaluate();
            bs.ulEvaluate();
            bothVariedSolution = [bothVariedSolution bs];
        end
    end

    mArray = {};
    for k = 1:2
        array = zeros(Problem.Du+Problem.Dl,Problem.Du + Problem.Dl);
        for i = 1:Problem.Du
            delta1 = upperVariedSolution(i).ulObj(k) - basicsolution.ulObj(k);
            for j = 1:Problem.Dl
                delta2 = bothVariedSolution((i-1)*Problem.Dl +j).ulObj(k) - lowerVariedSolution(j).ulObj(k);

                value = abs(delta1-delta2);
                array(i,j+Problem.Du) = value;
                array(Problem.Du + j,i) = value;
            end
        end
        mArray(k) = {array};
    end
    array1 = mArray{1} +mArray{2};
    
    %upper constraint
    upperCons = size(basicsolution.ulCon,2) - size(basicsolution.llCon,2);
    for t = 1:upperCons
        array = zeros(Problem.Du+Problem.Dl,Problem.Du + Problem.Dl);
        for i = 1:Problem.Du
            delta1 = upperVariedSolution(i).ulCon(t) - basicsolution.ulCon(t);
            if abs(delta1) < epsilon
                continue;
            end
            for j = 1:Problem.Dl
                delta2 = lowerVariedSolution(j).ulCon(t) - basicsolution.ulCon(t);
                value = abs(delta2);
                array(i,j+Problem.Du) = value;
                array(Problem.Du + j,i) = value;
            end
        end
        array1  = array1 + array;
    end

    group1 = [];
    for i = 1:Problem.Du
        for j = 1:Problem.Dl
            if array1(i,Problem.Du+j) >= epsilon && ~ismember(j,group1)
                group1 = [group1 j];
            end
        end
    end 

    ulRelate = {};
    for i = 1:Problem.Dl
        tmpRelate = [];
        for j = 1:Problem.Du
            if array1(i+Problem.Du,j) >= epsilon
                tmpRelate = [tmpRelate j];
            end
        end
        ulRelate{end+1} = tmpRelate;
    end 

    for k = 1:2
        array = zeros(Problem.Du,Problem.Dl);
        for i = 1:Problem.Du
            delta1 = upperVariedSolution(i).llObj(k) - basicsolution.llObj(k);
            
            for j = 1:Problem.Dl
                delta2 = bothVariedSolution((i-1)*Problem.Dl +j).llObj(k) - lowerVariedSolution(j).llObj(k);

                value = abs(delta1-delta2);
                array(i,j+Problem.Du) = value;
                array(Problem.Du+j,i) = value;
            end
        end
        mArray(k) = {array};
    end
    array2 = mArray{1} + mArray{2};

    %lower constraint
    lowerCons = size(basicsolution.llCon,2);
    for t = 1:lowerCons
        array = zeros(Problem.Du+Problem.Dl,Problem.Du + Problem.Dl);
        for i = 1:Problem.Du
            delta1 = upperVariedSolution(i).llCon(t) - basicsolution.llCon(t);
            if abs(delta1) < epsilon
                continue;
            end
            for j = 1:Problem.Dl
                delta2 = lowerVariedSolution(j).llCon(t) - basicsolution.llCon(t);
                value = abs(delta2);
                array(i,j+Problem.Du) = value;
                array(Problem.Du + j,i) = value;
            end
        end
        array2  = array2 + array;
    end

    group2 = [];
    for i = 1:Problem.Du
        for j = 1:Problem.Dl
            if array2(i,Problem.Du+j) >= epsilon && ~ismember(j,group2)
                group2 = [group2 j];
            end
        end
    end

    llRelate = {};
    for i = 1:Problem.Dl
        tmpRelate = [];
        for j = 1:Problem.Du
            if array2(i+Problem.Du,j) >= epsilon
                tmpRelate = [tmpRelate j];
            end
        end
        llRelate{end+1} = tmpRelate;
    end

    Relate = {};
    for i = 1:Problem.Dl
        if size(union(ulRelate{i},llRelate{i}),1) == 0
            Relate{end+1} = [];
        else
            Relate{end+1} = union(ulRelate{i},llRelate{i});
        end
    end

    archive4 = intersect(group1,group2);
    archive4 = getReal(archive4);

    archive3 = setdiff(group2,group1);
    archive3 = getReal(archive3);

    archive2 = setdiff(group1,group2);
    archive2 = getReal(archive2);

    archive1 = setdiff(setdiff([1:Problem.Dl],group1),group2);
    archive1 = getReal(archive1);

end

function archive = getReal(archive)
    if size(archive,1) == 0 || size(archive,2) == 0
        archive = [];
    end
end
