function [ values ] = makeSubstacks( stack, pixel_coordinates )
%% This function extracts x-y-z substacks from a TIFF file based on pixel-wise ROI coordinates.
% ----------------------------------------------------------------------
% Author: Dr. Knut Kirmse, 2018-2022.
% ----------------------------------------------------------------------
% INPUT:
% • stack: 3D matrix of a x-y-z TIFF stack
% • pixel_coordinates: 2D matrix generated in Fiji/ImageJ (macro: ReturnCoordinatesOfImageROI.ijm)
% containing coordinates (in MATLAB notation) of all pixels of each ROI for
% which a substack shall be extracted ...
% ... having the following structure:
%        col1       col2        col3        col4        ...     colN-1      colN
% row1   header     header      header      header      ...     header      header
% row2   ROI1,y1    ROI1,x1     ROI2,y1     ROI2,x1     ...     ROIn/2,y1   ROIn/2,x1
% row3   ROI1,y2    ROI1,x2     ROI2,y2     ROI2,x2     ...     ROIn/2,y2   ROIn/2,x2
% ...
% rowN   ROI1,yn    ROI1,xn     ROI2,yn     ROI2,xn     ...     ROIn/2,yn   ROIn/2,xn
% (Empty cells are filled with NaNs.)
% ----------------------------------------------------------------------
% OUTPUT:
% • values ... a cell array (dimensions: 2 x number of ROIs)
%       – row 1: ROI number
%       - row 2: substacks per ROI
% ----------------------------------------------------------------------
% Use:
% >> A = makeSubstacks(data, coord);
% ----------------------------------------------------------------------
%% Generate substacks
nROIs = 0.5 * size(pixel_coordinates, 2);
Nframes = size(stack, 3);
values = cell(2, nROIs);
for k = 1:nROIs
    values{1, k} = ['ROI_', int2str(k)];
    min_row = min(pixel_coordinates(:, k*2-1)); % minimum row of bounding rectangle
    max_row = max(pixel_coordinates(:, k*2-1)); % maximum row of bounding rectangle
    min_col = min(pixel_coordinates(:, k*2)); % minimum column of bounding rectangle
    max_col = max(pixel_coordinates(:, k*2)); % maximum column of bounding rectangle
    values{2, k} = stack(min_row:max_row, min_col:max_col, :);
    for m = min_row:max_row
        for n = min_col:max_col
            aux1 = pixel_coordinates(:, k*2-1) == m;
            aux2 = pixel_coordinates(:, k*2) == n;
            if sum(aux1 .* aux2) ~= 1
                values{2, k}(m-min_row+1, n-min_col+1, :) = NaN;
            end
        end
    end
end