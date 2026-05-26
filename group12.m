% Uploading data - .mat file?
data = load('modal_extraction.mat');

% Assembling variables
s = data.sensibility; 
raw_data = data.measures/s;  % Displacements, by a vibrometer i think

% Do we have 1 measuring point ??
%  we can repeat measurements more than 1 time

min_d = 0; % Looking at the data we chose a threshold, to discard
% the initial time of nothing

fs = data.fs;
dt = 1/fs;

for ii=1:length(raw_data(:,1)) % ! if they are column is 1,:
    index = find(d(ii,:) > min_d,'1','first');
    d(ii).data = d.raw_data(ii, index:end); % each data can have different 
    % dimension
    d(ii).t = 0:dt:dt*length(d(ii).data);
end

% We now have n repitition of the random signal we gave
% Problem: if the signal is random, we can't average to take out the noise
