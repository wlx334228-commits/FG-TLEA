classdef DS1 < PROBLEM
% <multi> <real> <large/none>
% Nu --- 20 ---ulPopulation size
% Nl --- 20 ---llPopulation size
% Mu --- 2 ---ulNumber of objectives
% Ml --- 2 ---llNumber of objectives
% Du --- 10 ---Number of decision variables
% Dl --- 10 ---Number of decision variables
% maxFE --- 1500000 --- Maximum number of function evaluations
    methods
        %% Default settings of the problem
        function Setting(obj)
            k= 10;
            ulDimMin = -k*ones(1,k);            
            ulDimMax = k*ones(1,k);
            ulDimMin(1) = max(ulDimMin(1),1);
            ulDimMax(1) = min(ulDimMax(1),4);
            obj.ulMin    = ulDimMin;
            obj.ulMax    = ulDimMax;
            
            obj.llMin    = -k*ones(1,k);
            obj.llMax    = k*ones(1,k); 
            
        end
        
        %% Calculate objective values
        function PopObj = ulCalObj(obj,ulPopDec,llPopDec)
            k= 10;
            r = 0.1;
            alpha = 1;
            gamma = 1;
            tao = 1;
            j = 2:k;
            PopObj(:,1) = (1 + r - cos(alpha * pi * ulPopDec(:,1))) + sum((ulPopDec(:,(j)) - (j-1)/2).^2) + tao * sum((llPopDec(:,(j)) - ulPopDec(:,(j))).^2)-...
                r * cos(gamma * (pi/2) * (llPopDec(:,1)/ulPopDec(:,1)));
            PopObj(:,2) = (1 + r - sin(alpha * pi * ulPopDec(:,1))) + sum((ulPopDec(:,(j)) - (j-1)/2).^2) + tao * sum((llPopDec(:,(j)) - ulPopDec(:,(j))).^2)-...
                r * sin(gamma * (pi/2) * (llPopDec(:,1)/ulPopDec(:,1)));
        end
        function PopObj = llCalObj(obj,ulPopDec,llPopDec)
            k= 10;
            PopObj(:,1) = llPopDec(:,1).^2 +sum((llPopDec(:,(2:k)) - ulPopDec(:,(2:k))).^2) + sum(10*(1-cos(pi/k*(llPopDec(:,(2:k))-ulPopDec(:,(2:k))))));
            PopObj(:,2) = sum((llPopDec - ulPopDec).^2) + sum(10 * abs(sin(pi/k*(llPopDec(:,(2:k))-ulPopDec(:,(2:k))))));
        end
    
        %% Sample reference points on Pareto front
        function P = GetOptimum(obj,N)
            r=0.1;
            P(:,1) = (0:(1+r)/(N-1):1+r)';
            P(:,2) =  -((1+r).^2-(P(:,1)-(1+r)).^2).^0.5 + (1+r);
        end
        %% Generate the image of Pareto front
        function R = GetPF(obj)
            R = obj.GetOptimum(500);
        end
        function R = GetName(obj)
            R = 'DS1';
        end
    end
end