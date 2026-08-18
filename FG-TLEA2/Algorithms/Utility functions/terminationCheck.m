function nofinish = terminationCheck(Archive,G,t,form)
    if nargin<4
        form =1;
    end
    if nargin<3
        t =500;
    end
    if isempty(Archive)
        nofinish = 1;
    else
        tau = size(Archive,2);
        mergePop = [];
        for i = 1:tau
            mergePop = [mergePop,Archive{1,i}];
        end
        if form ==1
            Objs = mergePop.llObjs;
        else
            Objs = mergePop.ulObjs;
        end
        nadirPoint = zeros(1,size(Objs,2));
        for i = 1:size(Objs,2)
            nadirPoint(1,i) = max(Objs(:,i));
        end
        hv = zeros(1,tau);
        for i = 1:tau
            if form ==1
                PopObj = Archive{1,i}.llObjs;
            else
                PopObj = Archive{1,i}.ulObjs;
            end
            hv(i) = llHV(PopObj,nadirPoint);
        end
        maxHV = max(hv);
        minHV = min(hv);
        if nargin<2 
            nofinish = ((maxHV - minHV) / (maxHV + minHV)>= 0.001);
        else
            nofinish = ((maxHV - minHV) / (maxHV + minHV)>= 0.001) * (G<t); %epsilonL
        end
    end
end
    
  
function score = llHV(PopObj,nadirPoint)
    offset = 0.1;
    RefPoint = nadirPoint+offset;
    N  = size(PopObj,1);
    if isempty(PopObj)
        score = 0;
    else
        pl = sortrows(PopObj);
        score = abs((RefPoint(1)-pl(1,1))*(RefPoint(2)-pl(1,2))); 
        for i = 2:N
            score = score + abs((RefPoint(1)-pl(i,1))*(pl(i-1,2)-pl(i,2)));
        end    
    end
end 

    
    