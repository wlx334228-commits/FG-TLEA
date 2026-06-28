classdef DS2 < PROBLEM
% <multi> <real> <large/none>
% Nu --- 20 ---ulPopulation size
% Nl --- 20 ---llPopulation size
% Mu --- 2 ---ulNumber of objectives
% Ml --- 2 ---llNumber of objectives
% Du --- 10 ---Number of decision variables
% Dl --- 10 ---Number of decision variables
% maxFE --- 3000000 --- Maximum number of function evaluations
    methods
        %% Default settings of the problem
        function Setting(obj)
            k= 10;
            ulDimMin = -k*ones(1,k);            
            ulDimMax = k*ones(1,k);
            ulDimMin(1) = max(ulDimMin(1),0.001);
            obj.ulMin    = ulDimMin;
            obj.ulMax    = ulDimMax;
            
            obj.llMin    = -k*ones(1,k);
            obj.llMax    = k*ones(1,k); 
            
        end
        
        %% Calculate objective values
        function PopObj = ulCalObj(obj,ulPopDec,llPopDec)
            k= 10;
            if ulPopDec(:,1)>=0 && ulPopDec(:,1)<=1
                v1 = cos(0.2*pi) * ulPopDec(:,1) + sin(0.2*pi) * sqrt(abs(0.02*sin(5*pi*ulPopDec(:,1))));
            elseif ulPopDec(:,1)>1
                v1 = ulPopDec(:,1) - (1 - cos(0.2*pi));
            end

            if ulPopDec(:,1)>=0 && ulPopDec(:,1)<=1
                v2 = -sin(0.2*pi) * ulPopDec(:,1) + cos(0.2*pi) * sqrt(abs(0.02*sin(5*pi*ulPopDec(:,1))));
            elseif ulPopDec(:,1)>1
                v2 = 0.1*(ulPopDec(:,1) -1) - sin(0.2*pi);
            end

            r = 0.25;
            gamma = 4;
            tao = 1;
            PopObj(:,1) = v1 + sum(ulPopDec(:,2:k).^2 + 10 * (1 - cos(pi/k * ulPopDec(:,2:k)))) + tao *...
    sum((llPopDec(:,2:k) - ulPopDec(:,2:k)).^2) - r * cos(gamma* (pi/2) * (llPopDec(:,1)/ulPopDec(:,1)));
            PopObj(:,2) = v2 + sum(ulPopDec(:,2:k).^2 + 10 * (1 - cos(1 - cos(pi/k * ulPopDec(:,2:k))))) + tao *...
    sum((llPopDec(:,2:k) - ulPopDec(:,2:k)).^2) - r* sin(gamma *(pi/2)* (llPopDec(:,1)/ulPopDec(:,1)));
        end
        function PopObj = llCalObj(obj,ulPopDec,llPopDec)
            k= 10;
            PopObj(:,1) = llPopDec(:,1).^2 + sum((llPopDec(:,2:k) - ulPopDec(:,2:k)).^2);
            PopObj(:,2) = sum((1:k).*(llPopDec - ulPopDec).^2);
        end
        %% Sample reference points on Pareto front
        function P = GetOptimum(obj,N)
            r=0.25;
            y1=[0.001,0.2,0.4,0.6,0.8,1];
            u1=cos(0.2*pi)*y1+sin(0.2*pi)*sqrt(abs(0.02*sin(5*pi*y1)));
            u2=-sin(0.2*pi)*y1+cos(0.2*pi)*sqrt(abs(0.02*sin(5*pi*y1)));
            PF1=[];PF2=[];
            theta=0:2*pi/(N-1):2*pi;
            for i =1:length(y1)
                Circle1=u1(i)+r*cos(theta);
                Circle2=u2(i)+r*sin(theta);
                PF1=[PF1,Circle1];
                PF2=[PF2,Circle2];
            end
            pf=([PF1;PF2])';
            NonDominated = NDSort(pf,1) == 1;
            P=pf(NonDominated',:);
        end
        %% Generate the image of Pareto front
        function R = GetPF(obj)
            R = obj.GetOptimum(565);
        end
        function R = GetName(obj)
            R = 'DS2';
        end
    end
end