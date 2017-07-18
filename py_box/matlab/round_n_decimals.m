function out = round_n_decimals(in, n, varargin)
    %Parse inputs
    p = inputParser;
    defaultRoundType = 'round';
    addRequired(p, 'in');
    addOptional(p, 'n', 0);
    addOptional(p, 'roundType', defaultRoundType);
    parse(p, in, n, varargin{:});

    %Round the value
    switch roundType
        case 'round'
            out = round(in*10^n)/10^n;
        case 'ceil'
            out = ceil(in*10^n)/10^n;
        case 'floor'
            out = floor(in*10^n)/10^n;
        otherwise
            error('Invalid roundType: %s.', roundType);
    end
end