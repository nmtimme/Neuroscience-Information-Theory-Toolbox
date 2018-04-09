%% INVERSEFNOISE - generates 1/f^alpha noise
% Generates noise with 1/f^alpha spectral power density over a specified
% range of frequencies. 
%
% Syntax: [noise] = inverseFNoise(n, alpha, Fs, varargin)
%
% Input:
%   n (integer) - the number of desired data points in the sample.
%   alpha (double) - the exponent of the noise spectrum (1/(f^alpha)).
%   Fs (double) - sampling frequency in Hz. 
%
% Optional Inputs:
%   
%   inverseFNoise(..., 'minF', minF) (double) - the minimum frequency in 
%       Hz. (default = the  minimum  frequency possible given n and Fs 
%       (Fs/n).)
%   inverseFNoise(..., 'maxF', maxF) (double) - the maximum frequency in 
%       Hz. (default = the maximum frequency possible given n and Fs 
%       (Fs/2).)
%   inverseFNoise(..., 'constLow' , constLow) (double) - defines the
%       frequency range below minF to include a constant power. The range
%       is defined in powers of 10. (default = 0)
%
% Outputs:
%   noise (1 by n double array) - the 1/f^alpha noise
%
% Examples:
%   10 seconds of 1/f noise with no cut offs sampled at 1 kHz
%       [noise] = inverseFNoise(10000, 1, 1000);
%   10 seconds of 1/f noise sampled at 1 kHz with no noise below 10 Hz.
%       [noise] = inverseFNoise(10000, 1, 1000, 'minF', 10);
%   100 seconds of 1/(f^1.5) noise sampled at 10 kHz with no noise below 1
%   Hz and constant noise power from 1 to 10 Hz.
%       [noise] = inverseFNoise(10^6, 1.5, 10000, 'minF', 1, 'constLow', 10);
%       
%   
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%

% Author: Nicholas Timme
% Email: nicholas.m.timme@gmail.com
% February 2016; Last revision: 5-Feb-2016


function [noise] = inverseFNoise(n, alpha, Fs, varargin)
%% Parse command line for parameters

% Calculate the minimum and maximum frequencies
minF = Fs/n;
maxF = Fs/2;
constLow = 0;

iVarArg = 1;
while iVarArg <= length(varargin)
    argOkay = true;
    switch varargin{iVarArg},
        case 'minF',        minF = varargin{iVarArg+1}; iVarArg = iVarArg + 1;
        case 'maxF',        maxF = varargin{iVarArg+1}; iVarArg = iVarArg + 1;
        case 'constLow',    constLow = varargin{iVarArg+1}; iVarArg = iVarArg + 1;
        otherwise,
            argOkay = false;
    end
    if ~argOkay
        disp(['(inverseFNoise) Ignoring invalid argument #' num2str(iVarArg+1)]);
    end
    iVarArg = iVarArg + 1;
end

%% Perform Initial Calculations

% Calculate the constant for the filter (assumes the maximum point on the
% filter is 1)
b = 10^(log10(1) + alpha*log10(minF));

% Calculate the frequencies
if rem(n,2) == 0
    bins = 1:(n/2);
    freqs = (Fs/n)*[0,bins,fliplr(bins(1:(end-1)))];
elseif rem(n,2) == 1
    bins = 1:floor(n/2);
    freqs = (Fs/n)*[0,bins,fliplr(bins)];
end

% Make the filter
filter = b*freqs.^(-alpha);
filter(freqs < 10^(log10(minF) - constLow)) = 0;
filter((freqs >= 10^(log10(minF) - constLow)) & (freqs < minF)) = 1;
filter(freqs > maxF) = 0;

%% Generate the Noise

% Generate white noise
noise = randn([1,n]);

% Fourier transform the noise
noise = fft(noise);

% Apply the 1/f^alpha filter
noise = filter.*noise;

% Inverse fourier transform the noise
noise = ifft(noise);

