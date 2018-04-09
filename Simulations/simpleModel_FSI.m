%% Fast Spiking Interneuron (FSI) Izhikevich Model Neuron
% Produces simulated data for a simple model of a fast spiking interneuron
% from Izhikevich 2007, page 299, Figure 8.27.
%
% Syntax: [mv,spkCt,spkLoc] = simpleModel_FSI(vInput,T,tau,PLOT,nzInput)
%
% Input:
%   vInput (1 by number of time steps double vector) - the current input to
%       the neuron in pA.
%   T (double) - the total time to run the model in ms.
%   tau (double) - the bin size in ms.
%   PLOT - if PLOT = 1, a plot of the voltage trace for the neuron will be
%       produced.
%   nzInput (1 by number of time steps double vector in units of mV) - 
%       membrane voltage noise.
%
% Outputs:
%   mV (double vector) - membrane potential of the model neuron in mV.
%   spkCt (integer) - number of spikes that occurred.
%   spkLoc (double vector) - a list of the time bins that contained spikes.
%
% Examples:
%   20 seconds of data binned at 1 ms with no noise and sinusoidal input
%   current
%       [mv,spkCt,spkLoc] = simpleModel_FSI(200*sin((1:(20000/0.1)).*(2*pi/10000)) + 200,20000,0.1,1,zeros([1,20000/0.1]));
%
% Citation: Izhikevich, E.M., 2007. Dynamical systems in neuroscience. MIT
% Press, Cambridge.
%   
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%

% Author: Nicholas Timme, Christopher Lapish
% Email: nicholas.m.timme@gmail.com
% February 2016; Last revision: 5-Dec-2016


%%FS - Layer 5 rat visual cortex fast-spiking (FS) interneuron (fig8.27 from 2007 book)
function [mv,spkCt,spkLoc] = simpleModel_FSI(vInput,T,tau,PLOT,nzInput)
%% FS interneuron parameters
C=20; vr=-55; vt=-40; k=3.5;
a=0.2; b=0.025; c=-45; vb=-55;
vpeak=25;

%% Prepare Variables

% Number of samples
n = round(T/tau);

% Voltage and recovery current
mv = vr*ones(1,n);
u = 0*mv;

% Spike records
spkCt = 0;
spkLoc = [];
m = 1;

%% Run the Simulation
for i=1:n-1                         % forward Euler method
    mv(i+1)=(mv(i)+tau*(k*(mv(i)-vr)*(mv(i)-vt)-u(i)+vInput(i))/C)+nzInput(i);
    % For FS neurons, include nonlinear U(v): U(v) = 0 when v<vb ; U(v) = a(v-vb) when v>=vb
    if (mv(i+1) < vb)
        u(i+1) = u(i) + tau*a*(0-u(i));
    else
        u(i+1) = u(i) + tau*a*((b*(mv(i)-vb).^3)-u(i));
    end
    %  Check if spike occurred and need to reset
    if mv(i+1)>=vpeak
        mv(i)=vpeak;
        mv(i+1)=c;
        spkCt=spkCt+1;     % count number of spikes and location
        spkLoc=[spkLoc m];
    end
    m=m+1;  
end

if PLOT==1;
    %% Plotting tools
    figure
    subplot(2,1,1);
    [ax,h1,h2] = plotyy(tau*(1:n),mv,tau*(1:n),vInput);
    xlabel('time (ms)');
    axes(ax(1)); ylabel('membrane potential (mV)');
    axes(ax(2)); ylabel('input current (pA)');
    subplot(2,1,2);
    plot(mv,u);
    xlabel('membrane potential (mV)');
    ylabel('recovery variable, u');
end;
  
