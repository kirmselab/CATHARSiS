function [fluor, spikes] = sim_fluor(firing_rate, frame_rate, rise_tau, decay_tau, Nframes, seed)
%% This function simulates a Poisson spike train and its corresponding fluorescence trace.
% ----------------------------------------------------------------------
% Author: Dr. Knut Kirmse, 2019-2022.
% Based on sim_fluor_UFARSA by Vahid Rahmati (December, 2017).
% ----------------------------------------------------------------------
% Input:
% firing_rate: in [Hz], expected mean firing rate of the simulated neuron
% frame_rate: in [Hz], sampling frequency
% rise_tau: in [sec], rise time-constant of calcium transients
% decay_tau: in [sec], decay time-constant of calcium transients
% nFrames: number of simulated frames
% seed: random seed; using a different seed will lead to a different spike train
% ----------------------------------------------------------------------
% Output:
% fluor: simulated fluorescence trace, acquired at the user-determined sampling frequency
% spikes: the underlying simulated Poisson spike-count at the given sampling frequency
% ----------------------------------------------------------------------
%% Compute
if ~isempty(seed)
    rng(seed)
end
si = 1 / frame_rate;  % sampling interval
spikes = poissrnd(firing_rate * si * ones(Nframes,1)); % spike train simulation
a = si:si:si*Nframes;
b = (1 - exp(-a ./ rise_tau)) .* exp(-a ./ decay_tau);
b = [zeros(1, numel(b)), b];
b = b / max(b); % template with amplitude = 1
fluor = conv(spikes, b, 'same'); % fluorescence trace
fluor = fluor';
spikes = spikes';
end