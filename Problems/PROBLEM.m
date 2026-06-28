classdef PROBLEM < handle & matlab.mixin.Heterogeneous
%PROBLEM - The superclass of problems.
%
%   This is the superclass of problems. An object of PROBLEM stores all the
%   settings of the problem.
%
% PROBLEM properties:
%   Nu               <read-only> population size--> upper level
%   Nl               <read-only> population size--> lower level
%   Mu               <read-only> number of objectives--> upper level
%   Ml               <read-only> number of objectives--> lower level
%   Du               <read-only> number of decision variables--> upper level
%   Dl               <read-only> number of decision variables--> lower level
%   maxFE           <read-only> maximum number of function evaluations--> upper+lower
%   ulFE              <read-only> number of consumed function evaluations--> upper level
%   llFE              <read-only> number of consumed function evaluations--> lower level
%   ulMin           <read-only> lower bound of decision variables--> upper level
%   ulMax           <read-only> upper bound of decision variables--> upper level
%   llMin           <read-only> lower bound of decision variables--> lower level
%   llMax           <read-only> upper bound of decision variables--> lower level
%   optimum         <read-only> optimal values of the problem
%   PF              <read-only> image of Pareto front--> upper level
%   gen             <read-only> current generation
%   encoding        <read-only>
%   parameter       <read-only> other parameters of the problem
%
% PROBLEM methods:
%   PROBLEM         <protected> the constructor, which sets all the properties specified by user
%   Setting         <public>    default settings of the problem
%   Initialization 	<public>    generate initial solutions
%   ulCalDec          <public>    repair invalid solutions
%   llCalDec          <public>    repair invalid solutions
%   ulCalObj          <public>    calculate the objective values of solutions
%   llCalObj          <public>    calculate the objective values of solutions
%   ulCalCon          <public>    calculate the constraint violations of solutions
%   llCalCon          <public>    calculate the constraint violations of solutions
%   GetOptimum      <public>    generate the optimums of the problem
%   GetPF          	<public>    generate the image of Pareto front
%   DrawDec         <public>    display a population in the decision space
%   DrawObj         <public>    display a population in the objective space
%   Current         <static>    get or set the current PROBLEM object
%   ParameterSet	<protected>	obtain the parameters of the problem

%------------------------------- Copyright --------------------------------
% Copyright (c) 2021 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    properties
        
        Mu  = 2;  Ml  = 2;              % Number of objectives        
        ulFE = 0; llFE = 0;  gen = 0;                % Number of consumed function evaluations
        encoding = 'real';
    end
    properties(SetAccess = protected)
        Nu;   Nl;            	% Population size
        Du;   Dl;                       % Number of decision variables
        maxFE;  % Maximum number of function evaluations
        ulMin; llMin;                	% Lower bound of decision variables
        ulMax; llMax;                 	% Upper bound of decision variables
        optimum;                    	% Optimal values of the problem
        PF;                          	% Image of Pareto front
        parameter = {};                	% Other parameters of the problem
    end
    methods(Access = protected)
        function obj = PROBLEM(varargin)
        %PROBLEM - The constructor of PROBLEM.
        %
        %   Problem = proName('Name',Value,'Name',Value,...) generates an
        %   object with the properties specified by the inputs. proName is
        %   PROBLEM or a subclass of PROBLEM.
        %
        %   If proName is PROBLEM, the properties objFcn and conFcn should
        %   be specified to define the objective and constraint functions.
        %   If proName is a subclass of PROBLEM, the properties objFcn and
        %   conFcn are useless since the objective and constraint functions
        %   are defined in the methods proName.CalObj and proName.CalCon.
        %
        %   The properties M, D, maxFE, encoding, lower, and upper will be
        %   further revised in the method proName.Setting.
        %
        %   Example:
        %       Problem = PROBLEM('objFcn',@(x)sum(x,2))
        %       Problem = DTLZ2('M',5,'D',10)

            isStr = find(cellfun(@ischar,varargin(1:end-1))&~cellfun(@isempty,varargin(2:end)));
            for i = isStr(ismember(varargin(isStr),{'Nu','Nl','Mu','Ml','Du','Dl','maxFE','parameter'}))
                obj.(varargin{i}) = varargin{i+1};
            end
            obj.Setting();
            obj.optimum = obj.GetOptimum(10000);
            obj.PF      = obj.GetPF();
        end
    end
    methods
        function Setting(obj)
        %Setting - Default settings of the problem.
        %
        %   This function is expected to be implemented in each subclass of
        %   PROBLEM, which will be called automatically.
        end
        function Population = Initialization(obj,Nu,Nl)
        %Initialization - Generate initial solutions.
        %
        %   P = obj.Initialization() randomly generates the decision
        %   variables of obj.N solutions and returns the SOLUTION objects.
        %
        %   P = obj.Initialization(N) generates N solutions.
        %
        %   Example:
        %       Population = Problem.Initialization()
        
            if nargin < 2
            	Nu = obj.Nu;
                Nl = obj.Nl;
            end
            

            Population = [];
            ulPopDec = unifrnd(repmat(obj.ulMin,Nu,1),repmat(obj.ulMax,Nu,1));
            for i = 1:Nu
                llPopDec = unifrnd(repmat(obj.llMin,Nl,1),repmat(obj.llMax,Nl,1));
                s = SOLUTION(ulPopDec(i,:),llPopDec,i);
                s = s.llEvaluate();
%                 % add
%                 s = llOptimizer(s);
                s = s.ulEvaluate();
                s = NDSelection(s,'l');
                Population= [Population,s];
            end 
        end
        
        function PopDec = ulCalDec(obj,PopDec)
        %CalDec - Repair invalid solutions.
        %
        %   Dec = obj.CalDec(Dec) repairs the invalid (not infeasible)
        %   decision variables of Dec.
        %
        %   An invalid solution indicates that it is out of the decision
        %   space, while an infeasible solution indicates that it does not
        %   satisfy all the constraints.
        %
        %   This function will be used when obj.decFcn is not specified.
        %
        %   Example:
        %       PopDec = Problem.CalDec(PopDec)

            DimMax = obj.ulMax;
            DimMin = obj.ulMin;
            PopDec = max(min(PopDec,repmat(DimMax,size(PopDec,1),1)),repmat(DimMin,size(PopDec,1),1));
           
        end
        function PopDec = llCalDec(obj,PopDec)
            DimMax = obj.llMax;
            DimMin = obj.llMin;
            PopDec = max(min(PopDec,repmat(DimMax,size(PopDec,1),1)),repmat(DimMin,size(PopDec,1),1));
           
        end
        
        function ulPopObj = ulCalObj(obj,ulPopDec,llPopDec)
        %CalObj - Calculate the objective values of solutions.
        %
        %   Obj = obj.CalObj(Dec) returns the objective values of Dec.
        %
        %   This function will be used when obj.objFcn is not specified.
        %
        %   Example:
        %       PopObj = Problem.CalObj(PopDec)

            ulPopObj = zeros(size(ulPopDec,1),obj.Mu);
        end
        function llPopObj = llCalObj(obj,ulPopDec,llPopDec)
        
            llPopObj = zeros(size(llPopDec,1),obj.Ml);
        end
        
        function ulPopCon = ulCalCon(obj,ulPopDec,llPopDec)
        %CalCon - Calculate the constraint violations of solutions.
        %
        %   Con = obj.CalCon(Dec) returns the constraint violations of Dec.
        %
        %   This function will be used when obj.conFcn is not specified.
        %
        %   Example:
        %       PopCon = Problem.CalCon(PopDec)
        
            ulPopCon = zeros(size(ulPopDec,1),1);
        end
        function llPopCon = llCalCon(obj,ulPopDec,llPopDec)
        
            llPopCon = zeros(size(llPopDec,1),1);
        end
        
        function R = GetOptimum(obj,N)
        %GetOptimum - Generate the optimums of the problem.
        %
        %   R = obj.GetOptimum(N) returns N optimums of the problem for
        %   metric calculation.
        %
        %   For multi-objective optimization problems, an optimum can be a
        %   point on the Pareto front; if the Pareto front is unknown, an
        %   optimum can be a reference point for hypervolume calculation.
        %
        %   For multi-modal multi-objective optimization problems, an
        %   optimum can be the decision variables of a Pareto optimal
        %   solution.
        %
        %   For single-objective optimization problems, an optimum can be
        %   the minimum objective value of the problem.
        %
        %   Example:
        %       R = Problem.GetOptimum(10000)
        
            if obj.Mu > 1
                R = ones(1,obj.Mu);
            else
                R = 0;
            end
        end
        function R = GetPF(obj)
        %GetPF - Generate the image of Pareto front.
        %
        %   R = obj.GetPF() returns the image of Pareto front for objective
        %   visualization.
        %
        %   For bi-objective optimization problems, the image should be a
        %   one-dimensional curve.
        %
        %   For tri-objective optimization problems, the image should be a
        %   two-dimensional surface.
        %
        %   For constrained multi-objective optimization problems, the
        %   image can be the feasible region.
        %
        %   Example:
        %       R = Problem.GetPF()
        
            R = [];
        end
        function DrawDec(obj,Population)
        %DrawDec - Display a population in the decision space.
        %
        %   obj.DrawDec(P) displays the decision variables of population P.
        %
        %   Example:
        %       Problem.DrawDec(Population)
        
            Draw(cat(2,Population.ulDecs,Population.llDecs),{'\it x\rm_1','\it x\rm_2','\it x\rm_3'});
            
        end
        function DrawObj(obj,Population)
        %DrawObj - Display a population in the objective space.
        %
        %   obj.DrawObj(P) displays the objective values of population P.
        %
        %   Example:
        %       Problem.DrawObj(Population)

            ax = Draw(Population.ulObjs,{'\it f\rm_1','\it f\rm_2','\it f\rm_3'});
            if ~isempty(obj.PF)
                if ~iscell(obj.PF)
                    if obj.Mu == 2
                        plot(ax,obj.PF(:,1),obj.PF(:,2),'-k','LineWidth',1);
                    elseif obj.Mu == 3
                        plot3(ax,obj.PF(:,1),obj.PF(:,2),obj.PF(:,3),'-k','LineWidth',1);
                    end
                else
                    if obj.Mu == 2
                        surf(ax,obj.PF{1},obj.PF{2},obj.PF{3},'EdgeColor','none','FaceColor',[.85 .85 .85]);
                    elseif obj.Mu == 3
                        surf(ax,obj.PF{1},obj.PF{2},obj.PF{3},'EdgeColor',[.8 .8 .8],'FaceColor','none');
                    end
                    set(ax,'Children',ax.Children(flip(1:end)));
                end
            elseif size(obj.optimum,1) > 1 && obj.Mu < 4
                if obj.Mu == 2
                    plot(ax,obj.optimum(:,1),obj.optimum(:,2),'.k');
                elseif obj.Mu == 3
                    plot3(ax,obj.optimum(:,1),obj.optimum(:,2),obj.optimum(:,3),'.k');
                end
            end
        end
    end
    methods(Static, Sealed)
        function obj = Current(obj)
        %Current - Get or set the current PROBLEM object.
        %
        %   Pro = PROBLEM.Current() returns the current PROBLEM object.
        %
        %   PROBLEM.Current(Pro) sets the current PROBLEM object to Pro and
        %   sets Pro.evaluated to 0.
        %
        %   Example:
        %       Problem = PROBLEM.Current()
        
            persistent Problem;
            if nargin > 0
                Problem = obj;
            end
            if nargout > 0
                obj = Problem;
            end
        end
    end
	methods(Access = protected, Sealed)
        function varargout = ParameterSet(obj,varargin)
        %ParameterSet - Obtain the parameters of the problem.
        %
        %   [p1,p2,...] = obj.ParameterSet(v1,v2,...) sets the values of
        %   parameters p1, p2, ..., where each parameter is set to the
        %   value given in obj.parameter if obj.parameter is specified, and
        %   set to the value given in v1, v2, ... otherwise.
        %
        %   Example:
        %       [p1,p2,p3] = obj.ParameterSet(1,2,3)

            varargout = varargin;
            specified = ~cellfun(@isempty,obj.parameter);
            varargout(specified) = obj.parameter(specified);
        end
    end
end