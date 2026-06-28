function value = FGTLEAGetParameter(name,defaultValue)
%FGTLEAGetParameter - Read an optional global experiment parameter.
%
%   Parameters are stored in the global struct FGTLEA_PARAMS so experiment
%   scripts can change algorithm settings without editing core source files.

    global FGTLEA_PARAMS;
    value = defaultValue;

    if isempty(FGTLEA_PARAMS) || ~isstruct(FGTLEA_PARAMS)
        return;
    end

    if isfield(FGTLEA_PARAMS,name) && ~isempty(FGTLEA_PARAMS.(name))
        value = FGTLEA_PARAMS.(name);
    end
end
