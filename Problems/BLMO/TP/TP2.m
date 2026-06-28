classdef TP2 < PROBLEM
% <multi> <real> <large/none>
% Nu --- 5 ---ulPopulation size
% Nl --- 60 ---llPopulation size
% Mu --- 2 ---ulNumber of objectives
% Ml --- 2 ---llNumber of objectives
% Du --- 1 ---Number of decision variables
% Dl --- 14 ---Number of decision variables
% maxFE --- 500000 --- Maximum number of function evaluations
    methods
        %% Default settings of the problem
        function Setting(obj)
            obj.ulMin    = -1*ones(1,obj.Du);
            obj.ulMax    = 2*ones(1,obj.Du);
            obj.llMin    = -1*ones(1,obj.Dl);
            obj.llMax    = 2*ones(1,obj.Dl);
        end
        
        %% Calculate objective values
        function PopObj = ulCalObj(obj,ulPopDec,llPopDec)
            PopObj(:,1) = (llPopDec(:,1)-1).^2 + sum(llPopDec(:,2:end).^2) + ulPopDec.^2;
            PopObj(:,2) = (llPopDec(:,1)-1).^2 + sum(llPopDec(:,2:end).^2) + (ulPopDec-1).^2 ;
        end
        function PopObj = llCalObj(obj,ulPopDec,llPopDec)
            PopObj(:,1) = llPopDec(:,1).^2 + sum(llPopDec(:,2:end).^2);
            PopObj(:,2) = (llPopDec(:,1)-ulPopDec).^2 + sum(llPopDec(:,2:end).^2);
        end
        
        %% Sample reference points on Pareto front
        function P = GetOptimum(obj,N)
            F2=0:0.5/(N-1):0.5;
            F12=F2/2+(1-sqrt(F2/2)).^2;  
            PF=[F12',F2'];
            P=PF;
        end
        function R = GetPF(obj)
            R = obj.GetOptimum(500);
        end
        function R = GetName(obj)
            R = 'TP2';
        end
    end
end