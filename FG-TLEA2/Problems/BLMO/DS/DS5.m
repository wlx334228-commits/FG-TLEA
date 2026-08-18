classdef DS5 < PROBLEM
% <multi> <real> <large/none>
% Nu --- 5 ---ulPopulation size
% Nl --- 40 ---llPopulation size
% Mu --- 2 ---ulNumber of objectives
% Ml --- 2 ---llNumber of objectives
% Du --- 1 ---Number of decision variables
% Dl --- 9 ---Number of decision variables
% maxFE --- 1500000 --- Maximum number of function evaluations
    methods
        %% Default settings of the problem
        function Setting(obj)
            k = 5;                        
            l = 4;
            obj.ulMin    = 1;
            obj.ulMax    = 2;
            
            obj.llMin    = [0,-(k+l)*ones(1,k+l-1)];
            obj.llMax    = [1,(k+l )*ones(1,k+l-1)]; 
            
        end
        
        %% Calculate objective values
        function PopObj = ulCalObj(obj,ulPopDec,llPopDec)
            k = 5;                      
            PopObj(:,1) = (1 - llPopDec(:,1))*(1 + sum(llPopDec(:,2:k).^2 )) * ulPopDec(:,1);
            PopObj(:,2) = llPopDec(:,1)*(1 + sum(llPopDec(:,2:k).^2 )) * ulPopDec(:,1);
        end
        function PopObj = llCalObj(obj,ulPopDec,llPopDec)
            k = 5;                        
            l = 4;
            PopObj(:,1) = (1 - llPopDec(:,1)) * (1 + sum(llPopDec(:,k+1:k+l).^2))* ulPopDec(:,1);
            PopObj(:,2) = llPopDec(:,1) * (1 + sum(llPopDec(:,k+1:k+l).^2)) * ulPopDec(:,1);
        end
        %% Calculate constraint violations
        function PopCon = ulCalCon(obj,ulPopDec,llPopDec)
            PopCon = (1 - llPopDec(:,1)) * ulPopDec(:,1) + 1/2*llPopDec(:,1)*ulPopDec(:,1) - 2 + 1/5 *ceil([5*(1 - llPopDec(:,1))* ulPopDec(:,1) + 0.2]);
            PopCon = -PopCon;
        end
        function PopCon = llCalCon(obj,ulPopDec,llPopDec)
            PopCon = (1 - llPopDec(:,1)) * ulPopDec(:,1) + 1/2*llPopDec(:,1)*ulPopDec(:,1) - 2 + 1/5 *ceil([5*(1 - llPopDec(:,1))* ulPopDec(:,1) + 0.2]);
            PopCon = -PopCon;
        end
        %% Sample reference points on Pareto front
        function P = GetOptimum(obj,N)
            y1=[1,1.2,1.4,1.6,1.8];
            x11=2*(1-1./y1);
            x12=2*(1-0.9./y1);
            x1=[x11;x12];
            pf1=[];
            pf2=[];
            for i=1:length(y1)
                X = (x1(1,i):1/(N-1):x1(2,i))';
                PF1=(1-X)*y1(i)';
                PF2=X*y1(i);
                pf1=[pf1;PF1];
                pf2=[pf2;PF2];
            end
            P=[pf1,pf2];
        end
        %% Generate the image of Pareto front
        function R = GetPF(obj)
            R = obj.GetOptimum(667);
        end
        function R = GetName(obj)
            R = 'DS5';
        end
    end
end