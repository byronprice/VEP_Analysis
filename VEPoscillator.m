function [times,x] = VEPoscillator(numImpulses,impulseSeparation,psi,omega0)
% VEPoscillator.m
%  Model the VEP as a harmonic oscillator with damping and a driving
%  impulse
%
%INPUT: numImpulses - number of impulses
%       impulseSeparation - time between impulses
%       psi - damping coefficient
%       omega0 - intrinsic frequency of the oscillator
%
%OUTPUT: times - timestamps in milliseconds
%        x - 'position' which in our case is microVolts
%
%Created: 2016/10/25
%  Byron Price
%Updated: 2016/10/25
% By: Byron Price

dt = 1/100; % time in ms
T = 1000; % total time in ms
% equation has form 
%  F(t) -kx -cdx/dt = m d2x/dt2
% rewritten as
%  d2x/dt2 +2(psi)(omega-not)*dx/dt+(omega-not)^2 = F(t)/m

% numerical solution
% psi = 0.1 - damping coefficient
% omega0 = 1/20; % frequency in 1/ms

times = 0:dt:T;

x = zeros(length(times),1);
xprime = zeros(length(times),1);

impulse = zeros(length(times),1);
impulseLen = 10; % make it last 10 milliseconds
magnitude = -10; % 100 mV magnitude
impulseDecay = 1;

startTime = 50;
for ii=1:numImpulses;
    count = 1;
    for tt=times
        if tt >= startTime && tt <= (startTime+impulseLen)
            impulse(count) = magnitude*exp(-impulseDecay*(tt-startTime));
        end
        count = count+1;
    end
    startTime = startTime+impulseSeparation;
end

x2prime = 0;
for ii=2:length(times)
    x2prime = impulse(ii-1)-2*psi*omega0*xprime(ii-1)-omega0.^2*x(ii-1);
    xprime(ii) = xprime(ii-1)+dt*x2prime;
    x(ii) = x(ii-1)+dt*xprime(ii);
end

figure();subplot(2,1,1);plot(times,x);title('Voltage vs. Time for VEP Oscillator');
xlabel('Time (milliseconds)');ylabel('Voltage (\muV)');legend(sprintf('Damping = %3.2f',psi));
subplot(2,1,2);plot(times,impulse);title('Driving Impulse vs. Time');
xlabel('Time (milliseconds)');ylabel('Voltage (\muV)');

end
