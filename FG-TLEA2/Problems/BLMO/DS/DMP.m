classdef DMP < PROBLEM
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
            llCount = 3;
            obj.ulMin    = zeros(1,ulCount);
            obj.ulMax    = 250*ones(1,ulCount);
            
            obj.llMin    = zeros(1,llCount);
            obj.llMax    = 70*ones(1,llCount); 
            
        end
        
        %% Calculate objective values
        function PopObj = ulCalObj(obj,ulPopDec,llPopDec)
            PopObj(:,1) = -(3.38*ulPopDec(:,1) + 7.78*ulPopDec(:,2) + 8.54*llPopDec(:,1) - 2.35*llPopDec(:,2) + 4.97*llPopDec(:,3));
            PopObj(:,2) = -(6.64*ulPopDec(:,1) + 4.26*ulPopDec(:,2) + 4.67*llPopDec(:,1) + 4.59*llPopDec(:,2) + 3.73*llPopDec(:,3));
        end
        function PopObj = llCalObj(obj,ulPopDec,llPopDec)
            PopObj(:,1) = -(4.47*ulPopDec(:,1) + 5.46*ulPopDec(:,2) - 6.23*llPopDec(:,1) - 4.78*llPopDec(:,2) + 7.34*llPopDec(:,3));
            PopObj(:,2) = -(5.34*ulPopDec(:,1) + 3.74*ulPopDec(:,2) + 9.45*llPopDec(:,1) + 6.37*llPopDec(:,2) + 5.45*llPopDec(:,3));
        end
        %% Calculate constraint violations
        function PopCon = ulCalCon(obj,ulPopDec,llPopDec)
            PopCon(:,1) = 4.55*ulPopDec(:,1) + 7.35*ulPopDec(:,2) + 9.65*llPopDec(:,1) - 6.23*llPopDec(:,2) + 4.24*llPopDec(:,3) - 987;
            PopCon(:,2) = -5.33*ulPopDec(:,1) - 1.35*ulPopDec(:,2) + 2.67*llPopDec(:,1) - 4.22*llPopDec(:,2) + 1.75*llPopDec(:,3) - 135;
            PopCon(:,3) = -2.11*ulPopDec(:,1) + 2.67*ulPopDec(:,2) + 4.34*llPopDec(:,1) + 9.26*llPopDec(:,2) + 8.33*llPopDec(:,3) - 830;
            PopCon(:,4) = 2.42*ulPopDec(:,1) + 7.43*ulPopDec(:,2) + 4.51*llPopDec(:,1) - 3.56*llPopDec(:,2) + 1.46*llPopDec(:,3) - 565;
        end
        function PopCon = llCalCon(obj,ulPopDec,llPopDec)
            PopCon(:,1) = 3.67*ulPopDec(:,1) - 7.84*ulPopDec(:,2) - 6.78*llPopDec(:,1) - 5.87*llPopDec(:,2) + 1.26*llPopDec(:,3) - 105;
            PopCon(:,2) = 4.34*ulPopDec(:,1) + 9.26*ulPopDec(:,2) + 8.33*llPopDec(:,1) - 2.11*llPopDec(:,2) - 2.67*llPopDec(:,3) - 830;
            PopCon(:,3) = 4.51*ulPopDec(:,1) - 3.56*ulPopDec(:,2) + 1.46*llPopDec(:,1) + 2.42*llPopDec(:,2) + 7.43*llPopDec(:,3) - 565;
            PopCon(:,4) = -(4.47*ulPopDec(:,1) + 5.46*ulPopDec(:,2) - 6.23*llPopDec(:,1) - 4.78*llPopDec(:,2) + 7.34*llPopDec(:,3));
            PopCon(:,5) = -(5.34*ulPopDec(:,1) + 3.74*ulPopDec(:,2) + 9.45*llPopDec(:,1) + 6.37*llPopDec(:,2) + 5.45*llPopDec(:,3));
            
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