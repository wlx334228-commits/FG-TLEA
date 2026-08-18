function Population = removeDupl(input)
    % Remove duplicates in the population
    Decs = cat(2,input.ulDecs,input.llDecs);
    [~,ia,~] = unique(Decs,'rows');
    Population = input(ia);
end