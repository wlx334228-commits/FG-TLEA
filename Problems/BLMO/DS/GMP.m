classdef GMP < PROBLEM
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
            ulCount = 2;
            llCount = 1;
            obj.ulMin    = [0 0];
            obj.ulMax    = [100 1];
            
            obj.llMin    = 0;
            obj.llMax    = 100; 
            
        end
        
        %% Calculate objective values
        function PopObj = ulCalObj(obj,ulPopDec,llPopDec)
            tau = ulPopDec(:,1);
            alpha_ = ulPopDec(:,2);
            q = llPopDec(:,1);

            PopObj(:,1) = -tau*q;
            PopObj(:,2) = q;
        end
        function PopObj = llCalObj(obj,ulPopDec,llPopDec)
            tau = ulPopDec(:,1);
            alpha_ = ulPopDec(:,2);
            q = llPopDec(:,1);
        
            PopObj(:,1) = -((100-q)*q - (q^2+q) - tau*q);
            PopObj(:,2) = q;
%             PopObj(:,1) = -(4.47*ulPopDec(:,1) + 5.46*ulPopDec(:,2) - 6.23*llPopDec(:,1) - 4.78*llPopDec(:,2) + 7.34*llPopDec(:,3));
%             PopObj(:,2) = -(5.34*ulPopDec(:,1) + 3.74*ulPopDec(:,2) + 9.45*llPopDec(:,1) + 6.37*llPopDec(:,2) + 5.45*llPopDec(:,3));
        end
        %% Calculate constraint violations
        function PopCon = ulCalCon(obj,ulPopDec,llPopDec)
            PopCon = [];
        end
        function PopCon = llCalCon(obj,ulPopDec,llPopDec)
            tau = ulPopDec(:,1);
            alpha_ = ulPopDec(:,2);
            q = llPopDec(:,1);
            PopCon(:,1) = -((100-q)*q - (q^2+q) - tau*q);
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