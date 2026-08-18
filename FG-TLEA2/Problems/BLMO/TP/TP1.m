classdef TP1 < PROBLEM
% <multi> <real> <large/none>
% Nu --- 5 ---ulPopulation size
% Nl --- 12 ---llPopulation size
% Mu --- 2 ---ulNumber of objectives
% Ml --- 2 ---llNumber of objectives
% Du --- 1 ---Number of decision variables
% Dl --- 2 ---Number of decision variables
% maxFE --- 500000 --- Maximum number of function evaluations
    methods
        %% Default settings of the problem
        function Setting(obj)
            obj.ulMin    = zeros(1,obj.Du);
            obj.ulMax    = ones(1,obj.Du);
            obj.llMin    = -ones(1,obj.Dl);
            obj.llMax    = ones(1,obj.Dl);
        end
        
        %% Calculate objective values
        function PopObj = ulCalObj(obj,ulPopDec,llPopDec)
            PopObj(:,1) = llPopDec(:,1)-ulPopDec;
            PopObj(:,2) = llPopDec(:,2);
        end
        function PopObj = llCalObj(obj,ulPopDec,llPopDec)
            PopObj(:,1) = llPopDec(:,1);
            PopObj(:,2) = llPopDec(:,2);
        end
        
        %% Calculate constraint violations
        function PopCon = ulCalCon(obj,ulPopDec,llPopDec)
            PopCon(:,1) = 1+llPopDec(:,1)+llPopDec(:,2);
            PopCon(:,2) = ulPopDec.^2-(llPopDec(:,1).^2+llPopDec(:,2).^2);
            PopCon = -PopCon;
        end
        function PopCon = llCalCon(obj,ulPopDec,llPopDec)
            PopCon = ulPopDec.^2-(llPopDec(:,1).^2+llPopDec(:,2).^2);
            PopCon = -PopCon;
        end
        
        %% Sample reference points on Pareto front
        function R = GetOptimum(obj,N)
            R(:,2) = (-1:1/(N-1):0)';
            t =sqrt(2*(R(:,2)+0.5).^2+0.5);
            R(:,1)=-1-R(:,2)-t;
        end
        
        function R = llGetOptimum(obj,ulPopDec,N)
            R(:,2) = (-ulPopDec:1/(N-1):0)';
            R(:,1) = -sqrt(ulPopDec.^2-R(:,2).^2);
        end
        %% Generate the image of Pareto front
        function R = GetPF(obj)
            R = obj.GetOptimum(500);
        end
        function R = GetName(obj)
            R = 'TP1';
        end
    end
end