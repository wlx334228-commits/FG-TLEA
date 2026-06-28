classdef SOLUTION < handle
%SOLUTION - The class of a solution.
%
%   This is the class of a solution. An object of SOLUTION stores all the
%   properties including decision variables, objective values, constraint
%   violations, and additional properties of a solution.
%
% SOLUTION properties:
%   ulDec         <public>     decision variables of the solution
%   llDec         <public>     decision variables of the solution
%   ulObj         <public>     objective values of the solution
%   llObj         <public>     objective values of the solution
%   ulCon         <public>     constraint violations of the solution
%   llCon         <public>     constraint violations of the solution
%   add           <public>     additional properties of the solution
%   NDu           <read-only>
%   NDl           <read-only>
%   CDu           <read-only>
%   CDl           <read-only>

% 
% SOLUTION methods:
%   SOLUTION	<public>        the constructor, which sets all the
%                               properties of the solution
%   ulEvaluate    <public>
%   llEvaluate    <public>
%   ulDecs        <public>      	get the matrix of decision variables of
%                               multiple solutions
%   llDecs        <public>      	get the matrix of decision variables of
%                               multiple solutions
%   ulObjs        <public>        get the matrix of objective values of
%                               multiple solutions
%   llObjs        <public>        get the matrix of objective values of
%                               multiple solutions
%   ulCons        <public>        get the matrix of constraint violations of
%                               multiple solutions
%   llCons        <public>        get the matrix of constraint violations of
%                               multiple solutions
%   adds        <public>        get the matrix of additional properties of
%                               multiple solutions
%   ulBest        <public>        get the feasible and nondominated solutions
%                               among multiple solutions
%   llBest        <public>        get the feasible and nondominated solutions
%                               among multiple solutions
%   NDus        <public> 
%   NDls        <public> 
%   CDus        <public> 
%   CDls        <public> 

%------------------------------- Copyright --------------------------------
% Copyright (c) 2021 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    properties%(SetAccess = private)
        ulDec; llDec       % Decision variables of the solution
    end
    properties  
        ulObj; llObj;        % Objective values of the solution
        ulCon; llCon;       % Constraint violations of the solution
        add;        % Additional properties of the solution
        NDu;NDl;CDu;CDl;
        title;
        rankl;
    end
    methods
        function obj = SOLUTION(ulPopDec,llPopDec,AddPro)
        %SOLUTION - The constructor of SOLUTION.
        %
        %   P = SOLUTION(Dec) creates an array of SOLUTION objects with
        %   decision variables of Dec, where the objective values and
        %   constraint violations are calculated automatically.
        %
        %   P = SOLUTION(Dec,AddPro) also sets the additional properties
        %   (e.g., velocity) of solutions to the values of AddPro.
        %
        %   Dec and AddPro are matrices, where each row denotes a solution
        %   and each column denotes a dimension of the decision variables
        %   or additional properties.
        %
        %   Example:
        %       Population = SOLUTION(PopDec)
            if nargin > 0
                if size(ulPopDec,1) ==1
                    ulPopDec = repmat(ulPopDec,size(llPopDec,1),1);
                end
                obj(1,size(llPopDec,1)) = SOLUTION;
                Problem = PROBLEM.Current();
                ulPopDec  = Problem.ulCalDec(ulPopDec);
                llPopDec  = Problem.llCalDec(llPopDec);
                for i = 1 : length(obj)
                    obj(i).ulDec = ulPopDec(i,:);
                    obj(i).llDec = llPopDec(i,:);
                end

                if nargin > 2
                    if size(AddPro,1) ==1
                        AddPro = repmat(AddPro,size(llPopDec,1),1);
                    end
                    for i = 1 : length(obj)
                        obj(i).add = AddPro(i,:);
                    end
                end
            end
        end
        
        function obj = ulEvaluate(obj)
            % obj.ulEvaluate calculate ul objective and constrain
            Problem = PROBLEM.Current();
            ulPopDec  = obj.ulDecs;
            llPopDec  = obj.llDecs;
            for i = 1 : length(obj)
                if isempty(obj(i).ulObj)
                    obj(i).ulObj = Problem.ulCalObj(ulPopDec(i,:),llPopDec(i,:));
                    obj(i).ulCon = Problem.ulCalCon(ulPopDec(i,:),llPopDec(i,:));
                    Problem.ulFE = Problem.ulFE + 1;
                end
            end
        end
        function obj = llEvaluate(obj)
            Problem = PROBLEM.Current();
            ulPopDec  = obj.ulDecs;
            llPopDec  = obj.llDecs;
            for i = 1 : length(obj)
                if isempty(obj(i).llObj)
                    obj(i).llObj = Problem.llCalObj(ulPopDec(i,:),llPopDec(i,:));
                    obj(i).llCon = Problem.llCalCon(ulPopDec(i,:),llPopDec(i,:));
                    Problem.llFE = Problem.llFE + 1;
                end
            end
        end
        
        function value = ulDecs(obj)
        %Decs - Get the matrix of decision variables of a population.
        %
        %   Dec = obj.Decs returns the matrix of decision variables of
        %   multiple solutions obj.
        
            value = cat(1,obj.ulDec);
        end
        function value = llDecs(obj)
            value = cat(1,obj.llDec);
        end
        
        function value = ulObjs(obj)
        %objs - Get the matrix of objective values of a population.
        %
        %   Obj = obj.objs returns the matrix of objective values of
        %   multiple solutions obj.
        
            value = cat(1,obj.ulObj);
        end
        function value = llObjs(obj)
            value = cat(1,obj.llObj);
        end
        
        function value = ulCons(obj)
        %cons - Get the matrix of constraint violations of a population.
        %
        %   Con = obj.cons returns the matrix of constraint violations of
        %   multiple solutions obj.
        
            value = cat(1,obj.ulCon);
        end
        function value = llCons(obj)
            value = cat(1,obj.llCon);
        end
        
        function value = adds(obj,AddPro)
        %adds - Get the matrix of additional properties of a population.
        %
        %   Add = obj.adds(AddPro) returns the matrix of additional
        %   properties of multiple solutions obj. If any solution in obj
        %   does not contain an additional property, it will be set to the
        %   default value specified in AddPro.
            if nargin > 1
                for i = 1 : length(obj)
                    if isempty(obj(i).add) || AddPro(i) ~= obj(i).add
                        obj(i).add = AddPro(i);
                    end
                end
            end
            value = cat(1,obj.add);
        end

        function value = titles(obj,TitlePro)
        %adds - Get the matrix of additional properties of a population.
        %
        %   Add = obj.adds(AddPro) returns the matrix of additional
        %   properties of multiple solutions obj. If any solution in obj
        %   does not contain an additional property, it will be set to the
        %   default value specified in AddPro.
            if nargin > 1
                for i = 1 : length(obj)
                    if isempty(obj(i).title) || TitlePro(i) ~= obj(i).title
                        obj(i).title = TitlePro(i);
                    end
                end
            end
            value = cat(1,obj.title);
        end

        function value = NDus(obj,NDu)
            if nargin > 1            
                for i = 1 : length(obj)
                    if isempty(obj(i).NDu) || NDu(i,:) ~= obj(i).NDu 
                        obj(i).NDu = NDu(i,:);
                    end
                end
            end
            value = cat(1,obj.NDu);
        end
        function value = NDls(obj,NDl)
            if nargin > 1 
                for i = 1 : length(obj)
                    if isempty(obj(i).NDl) || NDl(i,:) ~= obj(i).NDl 
                        obj(i).NDl = NDl(i,:);
                    end
                end
            end
            value = cat(1,obj.NDl);
        end
        function value = CDus(obj,CDu)
            if nargin > 1
                for i = 1 : length(obj)
                    if isempty(obj(i).CDu) || CDu(i,:) ~= obj(i).CDu 
                        obj(i).CDu = CDu(i,:);
                    end
                end
            end
            value = cat(1,obj.CDu);
        end
        function value = CDls(obj,CDl)
            if nargin > 1
                for i = 1 : length(obj)
                    if isempty(obj(i).CDl) || CDl(i,:) ~= obj(i).CDl 
                        obj(i).CDl = CDl(i,:);
                    end
                end
            end
            value = cat(1,obj.CDl);
        end

        function P = ulBest(obj)
        %best - Get the best solutions in a population.
        %
        %   P = obj.best returns the feasible and non-dominated solutions
        %   among multiple solutions obj. If the solutions have a single
        %   objective, the feasible solution with minimum objective value
        %   is returned.
        
            Feasible = find(all(obj.ulCons<=0,2));
            if isempty(Feasible)
                Best = [];
            elseif length(obj(1).ulObj) > 1
                Best = NDSort(obj(Feasible).ulObjs,1) == 1;
            else
                [~,Best] = min(obj(Feasible).ulObjs);
            end
            P = obj(Feasible(Best));
        end
        function P = llBest(obj)
            Feasible = find(all(obj.llCons<=0,2));
            if isempty(Feasible)
                Best = [];
            elseif length(obj(1).llObj) > 1
                Best = NDSort(obj(Feasible).llObjs,1) == 1;
            else
                [~,Best] = min(obj(Feasible).llObjs);
            end
            P = obj(Feasible(Best));
        end
    end
end