function sim_rings(dir, inR, outR, distrange)
%% This function generates simulated TIF movies for testing CATHARSiS.
% ----------------------------------------------------------------------
% Author: Dr. Knut Kirmse, 2019-2022.
% ----------------------------------------------------------------------
% Input:
% dir: output directory
% inR: inner radius of simulated cells [pixels]
% outR: outer radius of simulated cells [pixels]
% distrange: a scalar or vector specifying the offset(s) between cells [pixels]
% ----------------------------------------------------------------------
% Output:
% TIF file(s): image stack(s) (simulated Ca2+ imaging data)
% TXT file(s): containing ROI coordinates of simluated cell(s)
% ----------------------------------------------------------------------
%% Parameters
options_sim.noisescale = 1e4;
%
options_sim.firing_rate = 0.05;
% in [Hz], expected mean firing rate of the simulated neuron (lambda)
options_sim.frame_rate = 10;
% in [Hz], sampling frequency
options_sim.rise_tau = 0.1;
% in [sec], rise time-constant of calcium transients
options_sim.decay_tau = 1;
% in [sec], decay time-constant of calcium transients
options_sim.nFrames = 10000;
% number of simulated frames
%% Simulate cells
im = zeros(2*outR+1, max(distrange)+2*outR+1);
im1 = im;
im1(outR+1, outR+1) = 1;
im1 = bwdist(im1, 'euclidean') <= outR & bwdist(im1, 'euclidean') >= inR;
coord1 = find(im1); % ring only
im1_wholecell = im;
im1_wholecell(outR+1, outR+1) = 1;
im1_wholecell = bwdist(im1_wholecell, 'euclidean') <= outR;
coord1_wholecell = find(im1_wholecell); % ring & center
coord2 = NaN(length(coord1), numel(distrange));
coord2_wholecell = NaN(length(coord1_wholecell), numel(distrange));
im2 = zeros(2*outR+1, max(distrange)+2*outR+1, numel(distrange));
for k = 1:numel(distrange)
    im2(:, :, k) = im;
    im2(outR+1, outR+1+distrange(k), k) = 1;
    im2(:, :, k) = bwdist(im2(:, :, k), 'euclidean') <= outR & bwdist(im2(:, :, k), 'euclidean') >= inR;
    coord2(:, k) = find(im2(:, :, k)); % ring only
    im2_wholecell = im;
    im2_wholecell(outR+1, outR+1+distrange(k), k) = 1;
    im2_wholecell(:, :, k) = bwdist(im2_wholecell(:, :, k), 'euclidean') <= outR;
    coord2_wholecell(:, k) = find(im2_wholecell(:, :, k)); % ring & center
end
%% Simulate DF(t)
spikes = NaN(3, options_sim.nFrames); % two cells plus background
fluor = NaN(3, options_sim.nFrames);
for k = 1:3
    [fluor(k, :), spikes(k, :)] = sim_fluor(options_sim.firing_rate, options_sim.frame_rate, ...
        options_sim.rise_tau, options_sim.decay_tau, ...
        options_sim.nFrames, k); % <=== cave: seed
end
fluor(3, :) = fluor(3, :) .* 0.5; % scaling of background
%% Generate F(xyt), add Poisson noise and save TIF stacks and ROI TXT files for all pairs of simulated cells
simdata = single(NaN(size(im, 1), size(im, 2), options_sim.nFrames));
index = zeros(length(coord1_wholecell)+1, 4); % ROI read in runCATHARSiS.m starts at row 2
im_bg = ones(size(im, 1), size(im, 2)); % background (i.e. all pixels within frame)
bg_offset = 1; % baseline offset of background
for m = 1:numel(distrange)
    % Generate F(xyt)
    for k = 1:options_sim.nFrames
        simdata(:, :, k) = (im1 .* fluor(1, k)) + (im2(:, :, m) .* fluor(2, k)) + ...
            im1 + im2(:, :, m) + bg_offset + (im_bg .* fluor(3, k));
    end
    % Add Poisson noise
    simdata = single(imnoise(simdata / options_sim.noisescale, 'poisson') * options_sim.noisescale);
    % Save TIF stack file
    out = Fast_Tiff_Write([dir, '\simdata_', sprintf('%03d', distrange(m)), '.tif'], 1, 0);
        % arg2 - pixel size, arg3 - compression
    for k = 1:options_sim.nFrames
        aux = permute(simdata(:,:,k), [2, 1, 3]);
        out.WriteIMG(aux);
    end
    out.close
    % Save TXT ROI file (ring & center)
    [index(2:end, 1), index(2:end, 2)] = ind2sub([size(im, 1), size(im, 2)], coord1_wholecell);
    [index(2:end, 3), index(2:end, 4)] = ind2sub([size(im, 1), size(im, 2)], coord2_wholecell(:, m));
    dlmwrite([dir, '\simdata_ROI_', sprintf('%03d', distrange(m)), '.txt'], index, '\t');
end