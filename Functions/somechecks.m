

function somechecks(options)
% Some checks of the specified parameters
isinteger = @(x) mod(x, 1) == 0;

if options.TemplatesExist ~= 0 && options.TemplatesExist ~= 1
    error('options.TemplatesExist must be 0 or 1')
end

if options.firstframe ~= 0 && options.firstframe ~= 1
    error('options.firstframe must be 0 or 1')
end

if options.smooth ~= 0 && options.smooth ~= 1
    error('options.smooth must be 0 or 1')
end

if isinteger(options.Ntemplates) == 0 || options.Ntemplates < 1 || options.Ntemplates > 15
    error('options.Ntemplates must be a positive integer >= 1 AND <=15')
end

if isinteger(options.Ntemplates_min) == 0 || options.Ntemplates_min < 1
    error('options.Ntemplates_min must be a positive integer >= 1')
end

if options.Ntemplates_min > options.Ntemplates
    error('options.Ntemplates must larger than or equal to options.Ntemplates_min')
end

if isinteger(options.frameextension) == 0 || options.frameextension < 0
    error('options.frameextension must be an integer >= 0')
end

if isinteger(options.SG_length) == 0 || options.SG_length < 0
    error('options.SG_length must be an integer >= 0')
end

if isinteger(options.SG_order) == 0 || options.SG_order < 0
    error('options.SG_order must be an integer >= 0')
end

if options.SG_order >= options.SG_length
    error('options.SG_order must be <= options.SG_length - 1')
end

if options.smooth_CandCrit ~= 0 && options.smooth_CandCrit ~= 1
    error('options.smooth_CandCrit must be 0 or 1')
end

if isinteger(options.ROIexpansion) == 0 || options.ROIexpansion < 0
    error('options.Ntemplates must be a positive integer >= 0')
end

if options.Ntemplates_min == 1 && options.frameextension == 0
    error('Spatial TM: options.Ntemplates_min must be >= 2, if options.frameextension == 0')
end

end