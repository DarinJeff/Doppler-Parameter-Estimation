clear all;
close all;
clc;
%%
N = 1000;
fs = 1000;
L = 2;

P = 2*L;
Ts = 1/fs;
tau = 1/P * N * Ts;
t = Ts:Ts:N*Ts;
tt = 0:tau:N*Ts;
h = @rectpuls;
s = pulstran(t,tt,h,tau/10);
s = [s(N - fix(N/(P*20)):N) , s(1:N - fix(N/(P*20)) - 1)];
%-------------------------------------------------------------------------
%figure
subplot(3,2,1);
plot (t,s);
title('input pulse train');
ylabel('h(t - P*tau)');
xlabel('t (in sec)');
subplot(3,2,2);
plot(t,angle(s));
title('phase input');
%-------------------------------------------------------------------------
% attenuation, delay and doppler
% Note:      tl < N*Ts*9/(10*P) 
%       and Δtl > N/(10*P*fs)
al = [0.4 0.5];
tl = [0.08 0.13];
vl = [0.4  0.7];

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
%figure;
subplot(3,2,3);
plot(t,abs(x));
title('received signal');
ylabel('x(t)')
xlabel('t (in sec)')
subplot(3,2,4);
plot(t,angle(x));
title('phase received');
ylabel('angle s(t)')
xlabel('t (in sec)')
subplot(3,2,5);
plot(t,abs(xr));
title('approximated received signal');
ylabel('xr(t)');
xlabel('t (in sec)');
subplot(3,2,6);
plot(t,angle(xr));
title('approximated received signal phase');
ylabel('angle xr(t)');
xlabel('t (in sec)');
%-------------------------------------------------------------------------
     
xp = [xr(fix(N/P)+1:fix(2*N/P)) , zeros(1,N - fix(N/P))];
xw = fft(xp);
fx = (-length(xw)/2:length(xw)/2 - 1)*fs/length(xw);

hh = [ s(1:fix(N/P)), zeros(1,N - fix(N/P)) ];
hw = fft(hh);
%-------------------------------------------------------------------------
figure();
subplot(2,2,1);
plot(t,hh);
title('transmitted pulse');
ylabel('hh');
xlabel('t (in sec)');

subplot(2,2,2);
plot(fx/(2*pi),abs(fftshift(hw)));
title('fft of transmitted pulse');
ylabel('abs(hw)');
xlabel('w (in rad/s)');

subplot(2,2,3);
plot(t,abs(xp));
title('time shifted pth reflected signal');
ylabel('xp(t + P*tau)');
xlabel('t (in sec)');

subplot(2,2,4);
plot(fx/(2*pi),abs(fftshift(xw)));
title('fft of time shifted pth reflection');
ylabel('abs(xw)');
xlabel('w (in rad/s)');
%-------------------------------------------------------------------------
%%
xwk = zeros(1,2*L+1);
hwk = zeros(1,2*L+1);
max_ft = zeros(l, fix(2*pi/max(tl)));

for wo = 1:fix(2*pi/max(tl))
    good_wo = true;   
    zpw = zeros(1,2*L+1);
    
    for k = 1:2*L + 1
       if(hw(k*wo)==0)
           good_wo = false;     
           break;
       end
       xwk(k) = xw(k*wo);           % Equation 3.6
       hwk(k) = hw(k*wo);
       zpw(k) =  xwk(k)/hwk(k);     % Equation 3.7
    end
    
    if(~good_wo)
        wo
       continue;
    end
    
    F = zeros(L+1, L+1);
    for i = 1:L+1
        for j = 1:L+1
            F(i,j) = zpw(L+1+i-j);
        end
    end

    % annihilating filter  
    
    if(rank(F)<L | rank(F)==L+1) 
        continue;
    end
    
    A = null(F);
    zqr = roots(A);

    tq = angle(zqr);
    
    ft = zeros(1, L);
    for k = 1:2*L+1
        for l = 1:L
            ft(1, l) = ft(1, l) + zpw(k)*exp(-1j*k*wo*tq(l));
        end
    end
    max_ft(:,wo) = ft;
end

[~, wopt] = max(max_ft, [] ,2);

%% Calculating for wo = wopt1 and wo = wopt2
for w = 1:2
   for k = -4:4
     xwk(k+5) = xw((abs(k)*wopt(w))+1);         %% Equation 3.6
   end

  for k = -4:4
    hwk(k+5) = hw((abs(k)*wopt(w))+1);
  end

  for k = -4:4
    zpw(k+5) =  xwk(k+5)/hwk(k+5);     %% Equation 3.7
  end

an1 = ((zpw(8)*zpw(5))-(zpw(7)*zpw(6)))/((zpw(6)*zpw(6))-(zpw(5)*zpw(7)));        %% Equation 3.11
an2 = ((zpw(8)*zpw(6))-(zpw(7)*zpw(7)))/((zpw(5)*zpw(7))-(zpw(6)*zpw(6))); 

zqr = roots([1,an1,an2]);

tq(w) = abs(acos(real(zqr(w)))/(wopt(w)));
end

ft_opt1 = 0;
ft_opt2 = 0;
for k = 1:9
    
    ft_opt1 = ft_opt1 + (exp(-1j*(k-5)*wopt1*(tl(1))) * exp(1j*(k-5)*wopt1*tq(1)));  %% Equation 3.14
    ft_opt2 = ft_opt2 + (exp(-1j*(k-5)*wopt2*(tl(2))) * exp(1j*(k-5)*wopt2*tq(2)));
end


%% removing delay effect
dopp = zeros(1,25);
for i = 1:L
    rem = [x(1, (tl(i)*1000) + 1 : (tl(i)*1000) + 25 )];
    dopp = dopp + rem;
end

%dopp = upsample(dopp, 40);
dopp = [dopp , zeros(1,975)];
%dopp =repmat(dopp,1,40);
%plot (t,abs(dopp));

