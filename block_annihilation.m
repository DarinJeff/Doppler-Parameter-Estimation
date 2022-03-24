clear all;
close all;
clc;
% Change values in user input section, final delay estimates in Tq
%% User Input
% Number of samples, sampling frequency, targets, duty ratio of pulse
N = 4000;
fs = 2000;
L = 4;
D = 0.1;
%--------------------------------------------------------------------------
al = [0.4, 0.5, 0.7, 0.2];
tl = [0.05, 0.1, 0.125, 0.012];
vl = [0.04, 0.07, 0.09, 0.03];
%% Generating Transmitted and received pulses
%------------------------------------------
% to ensure that the divisions are proper
if(D*N/(4*L) ~= fix(D*N/(4*L)))
    error('Enter D such that N*D/(4*L) is an integer')
end
%------------------------------------------
P = 2*L;
Ts = 1/fs;
tau = 1/P * N * Ts;
t = Ts:Ts:N*Ts;
tt = 0:tau:N*Ts;
h = @rectpuls;
s = pulstran(t,tt,h,tau*D);
s = [s(N - fix(N*D/(P*2)):N) , s(1:N - fix(N*D/(P*2)) - 1)];

% received signal, 
x = zeros(1,N);
% countains 'L' delayed versions of s
srl = zeros(L,N);
% approximated recieved signal
xr = zeros(1,N);

for l = 1:L 
    srl(l,:) = [zeros(1,fix(tl(l)*fs)), s(1:N - fix(tl(l)*fs))];
    x = x + al(l)*srl(l,:).*exp(-1j*vl(l)*(t));           
end

for p = 1:P 
    for l = 1:L
        xr(fix((p-1)*N/P)+1:fix(p*N/P)) = xr(fix((p-1)*N/P)+1:fix(p*N/P)) + al(l)*srl(l,fix((p-1)*N/P)+1:fix(p*N/P))*exp(-1i*vl(l)*(p-1)*tau);
    end
end
%-------------------------------------------------------------------------
%figures;
subplot(3,2,1);
plot (t,s);
title('input pulse train');
ylabel('h(t - P*tau)');
xlabel('t (in sec)');
subplot(3,2,2);
plot(t,angle(s));
title('phase input');
subplot(3,2,3);
plot(t,abs(xr));
title('approximated received signal');
ylabel('xr(t)');
xlabel('t (in sec)');
subplot(3,2,4);
plot(t,angle(xr));
title('approximated received signal phase');
ylabel('angle xr(t)');
xlabel('t (in sec)');
%-------------------------------------------------------------------------  
% noise
sigmaN = 0.000;
xr = xr + randn(1,N)*sigmaN;
%-------------------------------------------------------------------------  
subplot(3,2,5);
plot(t,abs(xr));
title('approximated noisy received signal');
ylabel('xr(t)');
xlabel('t (in sec)');
subplot(3,2,6);
plot(t,angle(xr));
title('approximated noisy received signal phase');
ylabel('angle xr(t)');
xlabel('t (in sec)');
%-------------------------------------------------------------------------  
% Get the P reflected pulses individually
n = N/P;
xp = zeros(P,n);
xpw = zeros(P,n);
for p = 0:P-1
   xp(p+1,:) = xr(p*n + 1: (p+1)*n);   
   xpw(p+1,:) = fft(xp(p+1,:));
end

fxp = (-n/2:n/2 - 1)*fs/n;
tp = Ts:Ts:n*Ts;

% cropped transmitted pulse
hh = s(1:n);
hw = fft(hh);
%-------------------------------------------------------------------------
figure();
subplot(2,2,1);
plot(tp,hh);
title('transmitted pulse');
ylabel('hh');
xlabel('t (in sec)');

subplot(2,2,2);
plot(fxp/(2*pi),abs(fftshift(hw)));
title('fft of transmitted pulse');
ylabel('abs(hw)');
xlabel('w (in rad/s)');

subplot(2,2,3);
plot(tp,abs(xp(1,:)));
title('1st reflected signal');
ylabel('xp(t + P*tau)');
xlabel('t (in sec)');

subplot(2,2,4);
plot(fxp/(2*pi),abs(fftshift(xpw(1,:))));
title('fft of 1st reflection');
ylabel('abs(xw)');
xlabel('w (in rad/s)');
%% Delay Estimation by method of annihilation
K = 2*L + 1;
xpwk = zeros(P,K);
hwk = zeros(1,K);
Tq = zeros(min( fix(n*Ts/max(tl)), fix(n/K) ), l);

for wo = 1: min( fix(n*Ts/max(tl)), fix(n/K))
    good_wo = true;   
    zpwk = zeros(P,K);
    
    for k = 0:K-1
        if(hw(k*wo+1)==0)
            good_wo = false;     
            break;
        end
        hwk(k+1) = hw(k*wo+1);
        for p = 1:P
            xpwk(p,k+1) = xpw(p,k*wo+1);
            zpwk(p,k+1) =  xpwk(p,k+1)/hwk(k+1);
        end
    end
    
    if(~good_wo)
       disp('h[kwo] = 0');
       continue;
    end
    
    F = zeros((L+1)*P, L+1);
    for p = 1:P
        for i = 1:L+1
            for j = 1:L+1
                F(i + (p-1)*(L+1),j) = zpwk(p, L+1+i-j);
            end
        end
    end

    % annihilating filter  


    % eigenvector with smallest eigenvalue
    [A,D] = eigs(F'*F,1,'SM');
    zqr = roots(A);
    
    % Ensuring that the angle is from 0 to -2pi
    arg = 2*pi - mod(angle(zqr), 2*pi);
    tq = arg*n*Ts/(2*pi*wo);
    
    Tq(wo, :) = tq;
end
% Best estimate is in last non zero row - wo
row_is_zero = all(Tq~=0,2);
wo = find(row_is_zero,1,'last');
tq = Tq(wo, :)
