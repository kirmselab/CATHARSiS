function [ DETECTION_CRITERION ] = SpatialTemplateMatch( stack_ROI, template_ROI , driftindex)
%% This function performs SPATIAL TEMPLATE MATCHING.
% The algorithm computes a detection criterion as previously
% defined for the temporal domain (Clements & Bekkers, Biophys J, 1997).
% ---------------------------------------------------------------------------
% Author: Dr. Knut Kirmse, 2019-2022.
% ---------------------------------------------------------------------------
% INPUT:
% • stack_ROI: a stack with dimensions xyt
% • template_ROI: an image with dimensions xy
% • driftindex: a logical vector of length t
%   (0 - no drift, 1 - drift, i.e. exclude frame)
% ---------------------------------------------------------------------------
% OUTPUT:
% • DETECTION_CRITERION: a vector of length t (NaN for excluded frames)
% ---------------------------------------------------------------------------
%% Compute DETECTION_CRITERION
if size(stack_ROI, 1) ~= size(template_ROI, 1) || size(stack_ROI, 2) ~= size(template_ROI, 2)
    error('XY dimensions of stack and template do not match.')
end
if size(stack_ROI, 3) ~= length(driftindex)
    error('The number of elements in driftindex does not match the number of frames in stack_ROI.')
end
finalFrames = find(driftindex == 0);
Nrows = size(stack_ROI, 1);
Ncols = size(stack_ROI, 2);
DETECTION_CRITERION = NaN(length(driftindex), 1);
%
template_ROI = reshape(template_ROI, [Nrows * Ncols, 1]);
template_ROI(any(isnan(template_ROI), 2), :) = [];
SUM_template = sum(template_ROI);
SUM_template_squared = sum(template_ROI .* template_ROI);
N = length(template_ROI);
for m = finalFrames
    % Compute SCALE
    DATA = stack_ROI(:, :, m);
    DATA = reshape(DATA, [Nrows * Ncols, 1]);
    DATA(any(isnan(DATA), 2), :) = [];
    SUM_template_DATA = sum(template_ROI .* DATA);
    SUM_DATA = sum(DATA);
    SCALE = (SUM_template_DATA - SUM_template * SUM_DATA / N) / (SUM_template_squared - SUM_template * SUM_template / N);
    % Compute OFFSET
    OFFSET = (SUM_DATA - SCALE * SUM_template) / N;
    % Compute FITTED_TEMPLATE
    FITTED_TEMPLATE = template_ROI * SCALE + OFFSET;
    % Compute STANDARD_ERROR
    SSE = sum((DATA - FITTED_TEMPLATE).^2);
    STANDARD_ERROR = sqrt(SSE / (N-1));
    % Compute DETECTION_CRITERION
    DETECTION_CRITERION(m) = SCALE / STANDARD_ERROR;
end
end