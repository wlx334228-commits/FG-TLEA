classdef DS3 < PROBLEM
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
            obj.ulMin    = zeros(1,k);
            obj.ulMax    = k*ones(1,k);
            
            obj.llMin    = -k*ones(1,k);
            obj.llMax    = k*ones(1,k); 
            
        end
        
        %% Calculate objective values
        function PopObj = ulCalObj(obj,ulPopDec,llPopDec)
            k= 10;
            tao = 1;
            R = 0.1 + 0.15 * abs(sin(2*pi*(ulPopDec(:,1)-0.1)));
            PopObj(:,1) = ulPopDec(:,1) + sum((ulPopDec(:,3:k)-(3:k)/2).^2) + tao*sum((llPopDec(:,3:k) - ulPopDec(:,3:k)).^2) -...
    R * cos(4 * atan((ulPopDec(:,2) - llPopDec(:,2))/(ulPopDec(:,1) - llPopDec(:,1))));
            PopObj(:,2) = ulPopDec(:,2) + sum((ulPopDec(:,3:k)-(3:k)/2).^2) + tao*sum((llPopDec(:,3:k) - ulPopDec(:,3:k)).^2) -...
    R * sin(4 * atan((ulPopDec(:,2) - llPopDec(:,2))/(ulPopDec(:,1) - llPopDec(:,1)))) ;
        end
        function PopObj = llCalObj(obj,ulPopDec,llPopDec)
            k= 10;
            PopObj(:,1) = llPopDec(:,1) + sum((llPopDec(:,3:k) - ulPopDec(:,3:k)).^2);
            PopObj(:,2) = llPopDec(:,2) + sum((llPopDec(:,3:k) - ulPopDec(:,3:k)).^2);
        end
        %% Calculate constraint violations
        function PopCon = ulCalCon(obj,ulPopDec,llPopDec)
            r = 0.2;
            PopCon(:,1) = ulPopDec(:,2) - (1 - ulPopDec(:,1).^2);
            PopCon(:,2) = -(llPopDec(:,1) - ulPopDec(:,1)).^2 - (llPopDec(:,2) - ulPopDec(:,2)).^2 + r.^2;
            PopCon = -PopCon;
        end
        function PopCon = llCalCon(obj,ulPopDec,llPopDec)
            r = 0.2;
            PopCon = -(llPopDec(:,1) - ulPopDec(:,1)).^2 - (llPopDec(:,2) - ulPopDec(:,2)).^2 + r.^2;
            PopCon = -PopCon;
        end
        %% Sample reference points on Pareto front
        function P = GetOptimum(obj,N)
            y1 = (0:0.1:1.5)';
            R = 0.1 + 0.15*abs(sin(2*pi*(y1-0.1)));
            y2=1-y1.^2;
            y2(y2<0)=0;
            PF1=[];PF2=[];
            theta=0:2*pi/(N-1):2*pi;
            for i =1:length(y1)
                Circle1=y1(i)+R(i)*cos(theta);
                Circle2=y2(i)+R(i)*sin(theta);
                PF1=[PF1,Circle1];
                PF2=[PF2,Circle2];
            end
            pf=[PF1;PF2]';
            NonDominated = NDSort(pf,1) == 1;
            P=pf(NonDominated',:);
        end
        %% Generate the image of Pareto front
        function R = GetPF(obj)
            R = obj.GetOptimum(361);
        end
        function R = GetName(obj)
            R = 'DS3';
        end

    end
end