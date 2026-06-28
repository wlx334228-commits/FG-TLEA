function llPop = llOptimizerArchive(llPop,archive1,archive2)
    Problem = PROBLEM.Current();
    ulDec = llPop(1).ulDec;
    add = llPop(1).add;
    llPop = NDSelection(llPop,'l',Problem.Nl);
    elitePop = llPop(llPop.NDls == 1);
    G = 0;
    llArchive = [];
    tau = 10;
    while terminationCheck(llArchive,G)
        llArchive = cell(1,tau);
        for i= 1:tau
            llOff = [];
            for j = 1:Problem.Nl
                DecP1 = llPop(j).llDec;
                DecP2 = elitePop((unidrnd(size(elitePop,2)))).llDec;
                while 1
                    randP = unidrnd(Problem.Nl);
                    if randP ~= j
                        DecP3 = llPop(randP).llDec;
                        break;
                    end
                end
                llOffDec = OperatorDE(DecP1,DecP2,DecP3,'l');
                llOffDec([archive1,archive2]) = DecP1([archive1,archive2]);
                Off = SOLUTION(ulDec,llOffDec,add);
                Off.llEvaluate();
                llOff = [llOff,Off];
            end
            llPop = NDSelection([llPop,llOff],'l',Problem.Nl);
            elitePop = llPop(llPop.NDls==1);
            llArchive{1,i} = elitePop;
            G=G+1;
        end
    end
end




