%--------------- Preparation
clear all;
close all;
clc;

tic;
% Time Domain 0 to T
T = 1000;
fs = 1/T;
t = (1:T)/T;
freqs = 2*pi*(t-0.5-1/T)/(fs);

% center frequencies of components
f_1 = 150;               % Doppler 1 
f_2 = 280;              % Doppler 2
f_3 = 420;              % Doppler 3
md_1 = 25;               % mD 1
md_2 = 15;              % mD 2
md_3 = 30;              % mD 3
al = [1, 0.5, 1.3];
bl = [0.3, 0.1, 0.7];

targets = 3;
x = 0;
f_t = zeros(1,targets);
f_t = [f_1 f_2 f_3];
md_t =[md_1 md_2 md_3];
% v_1 has all doppler terms
% v_2 has all mD terms
% v_2 expression is exp(1j*doppler*t).*exp(1j*amp of mD*sin(mD*t))

%for i = 1:targets
    v_1 = al(1)*(exp(1j*2*pi*f_1*t))+ al(2)*(exp(1j*2*pi*f_2*t)) + al(3)*(exp(1j*2*pi*f_3*t));
    v_2 = bl(1)*exp(1j*2*pi*f_1*t).*exp(1j*(1.5)*sin(2*pi*(md_1)*t)) + bl(2)*exp(1j*2*pi*f_2*t).*exp(1j*(0.5)*sin(2*pi*(md_2)*t)) + bl(3)*exp(1j*2*pi*f_3*t).*exp(1j*(0.25)*sin(2*pi*(md_3)*t));
    %x = x + v_1 + v_2;
%end

v_3 = 1/16*(exp(1j*2*pi*f_3*t));

% for visualization purposes
fsub = {};
wsub = {};
fsub{1} = v_1;
fsub{2} = v_2;
fsub{3} = v_3;
wsub{1} = 2*pi*f_1;
wsub{2} = 2*pi*f_2;
wsub{3} = 2*pi*f_3;

% composite signal, including noise
f = v_1 + v_2;% + 0.1*randn(size(v_2));
%f = awgn(f,1);
f_hat = fftshift((fft(f)));

% some sample parameters for VMD
alpha = 2000;        % moderate bandwidth constraint
tau = 0;            % noise-tolerance (no strict fidelity enforcement)
K = 3;              % 3 modes
DC = 0;             % no DC part imposed
init = 1;           % initialize omegas uniformly
tol = 1e-7;






%--------------- Run actual VMD code

[u, u_hat, omega] = VMD(f, alpha, tau, K, DC, init,f_t, tol);










%--------------- Visualization

% For convenience here: Order omegas increasingly and reindex u/u_hat
[~, sortIndex] = sort(omega(end,:));
omega = omega(:,sortIndex);
u_hat = u_hat(:,sortIndex);
u = u(sortIndex,:);
linestyles = {'b', 'g', 'm', 'c', 'c', 'r', 'k'};

figure('Name', 'Composite input signal' );
plot(t,f, 'k');
set(gca, 'XLim', [0 1]);

for sub = 1:length(fsub)
    figure('Name', ['Input signal component ' num2str(sub)] );
    plot(t,fsub{sub}, 'k');
    set(gca, 'XLim', [0 1]);
end

figure('Name', 'Input signal spectrum' );
loglog(freqs(T/2+1:end), abs(f_hat(T/2+1:end)), 'k');
set(gca, 'XLim', [1 T/2]*pi*2, 'XGrid', 'on', 'YGrid', 'on', 'XMinorGrid', 'off', 'YMinorGrid', 'off');
ylims = get(gca, 'YLim');
hold on;
for sub = 1:length(wsub)
    loglog([wsub{sub} wsub{sub}], ylims, 'k--');
end
set(gca, 'YLim', ylims);

figure('Name', 'Evolution of center frequencies omega');
for k=1:K
    semilogx(2*pi/fs*omega(:,k), 1:size(omega,1), linestyles{k});
    hold on;
end
set(gca, 'YLim', [1,size(omega,1)]);
set(gca, 'XLim', [2*pi,0.5*2*pi/fs], 'XGrid', 'on', 'XMinorGrid', 'on');

figure('Name', 'Spectral decomposition');

loglog(freqs(T/2+1:end)/(2*pi), abs(f_hat(T/2+1:end)), 'k:');
set(gca, 'XLim', [1 T/2]*pi*2, 'XGrid', 'on', 'YGrid', 'on', 'XMinorGrid', 'off', 'YMinorGrid', 'off');
hold on;
for k = 1:K
    loglog(freqs(T/2+1:end)/(2*pi), abs(u_hat(T/2+1:end,k)), linestyles{k});
end
set(gca, 'YLim', ylims);
title('frequency spectrum (log-log scale)');
ylabel('amplitude');
xlabel('frequency (Hz)');


for k = 1:K
    figure('Name', ['Reconstructed mode ' num2str(K)]);
    plot(t,u(k,:), linestyles{k});   hold on;
    if ~isempty(fsub)
        plot(t, fsub{min(k,length(fsub))}, 'k:');
    end
    set(gca, 'XLim', [0 1]);
end

toc;
time = toc;
%% vmd 

function [u, u_hat, omega] = VMD(signal, alpha, tau, K, DC, init,f_t, tol)
% Variational Mode Decomposition
% Authors: Konstantin Dragomiretskiy and Dominique Zosso
% zosso@math.ucla.edu --- http://www.math.ucla.edu/~zosso
% Initial release 2013-12-12 (c) 2013
%
% Input and Parameters:
% ---------------------
% signal  - the time domain signal (1D) to be decomposed
% alpha   - the balancing parameter of the data-fidelity constraint
% tau     - time-step of the dual ascent ( pick 0 for noise-slack )
% K       - the number of modes to be recovered
% DC      - true if the first mode is put and kept at DC (0-freq)
% init    - 0 = all omegas start at 0
%                    1 = all omegas start uniformly distributed
%                    2 = all omegas initialized randomly
% tol     - tolerance of convergence criterion; typically around 1e-6
%
% Output:
% -------
% u       - the collection of decomposed modes
% u_hat   - spectra of the modes
% omega   - estimated mode center-frequencies
%
% When using this code, please do cite our paper:
% -----------------------------------------------
% K. Dragomiretskiy, D. Zosso, Variational Mode Decomposition, IEEE Trans.
% on Signal Processing (in press)
% please check here for update reference: 
%          http://dx.doi.org/10.1109/TSP.2013.2288675



%---------- Preparations

% Period and sampling frequency of input signal
save_T = length(signal);
fs = 1/save_T;

% extend the signal by mirroring
T = save_T;
f_mirror(1:T/2) = signal(T/2:-1:1);
f_mirror(T/2+1:3*T/2) = signal;
f_mirror(3*T/2+1:2*T) = signal(T:-1:T/2+1);
f = f_mirror;

% Time Domain 0 to T (of mirrored signal)
T = length(f);
t = (1:T)/T;

% Spectral Domain discretization
freqs = t-0.5-1/T;

% Maximum number of iterations (if not converged yet, then it won't anyway)
N = 500;

% For future generalizations: individual alpha for each mode
Alpha = alpha*ones(1,K);

% Construct and center f_hat
f_hat = fftshift((fft(f)));
f_hat_plus = f_hat;
f_hat_plus(1:T/2) = 0;

% matrix keeping track of every iterant // could be discarded for mem
u_hat_plus = zeros(N, length(freqs), K);

% Initialization of omega_k
omega_plus = zeros(N, K);
switch init
    case 1
        for i = 1:K
            %omega_plus(1,i) = (0.5/K)*(i-1);
            omega_plus(1,i) = (fs)*f_t(i)*(1+randn(1)/9);
        end
    case 2
        omega_plus(1,:) = sort(exp(log(fs) + (log(0.5)-log(fs))*rand(1,K)));
    otherwise
        omega_plus(1,:) = 0;
end

% omega_plus(1,1) = 20;
% omega_plus(1,2) = 40;
% omega_plus(1,3) = 30;

% if DC mode imposed, set its omega to 0
%if DC
%    omega_plus(1,1) = 0;
%end

% start with empty dual variables
lambda_hat = zeros(N, length(freqs));

% other inits
uDiff = tol+eps; % update step
n = 1; % loop counter
sum_uk = 0; % accumulator



% ----------- Main loop for iterative updates



tic;
while ( uDiff > tol &&  n < N ) % not converged and below iterations limit
    
    % update first mode accumulator
    k = 1;
    sum_uk = u_hat_plus(n,:,K) + sum_uk - u_hat_plus(n,:,1);
    
    % update spectrum of first mode through Wiener filter of residuals
    u_hat_plus(n+1,:,k) = (f_hat_plus - sum_uk - lambda_hat(n,:)/2)./(1+Alpha(1,k)*(freqs - omega_plus(n,k)).^2);
    
    % update first omega if not held at 0
    if ~DC
        omega_plus(n+1,k) = (freqs(T/2+1:T)*(abs(u_hat_plus(n+1, T/2+1:T, k)).^2)')/sum(abs(u_hat_plus(n+1,T/2+1:T,k)).^2);
    end
    
    % update of any other mode
    for k=2:K
        
        % accumulator
        sum_uk = u_hat_plus(n+1,:,k-1) + sum_uk - u_hat_plus(n,:,k);
        
        % mode spectrum
        u_hat_plus(n+1,:,k) = (f_hat_plus - sum_uk - lambda_hat(n,:)/2)./(1+Alpha(1,k)*(freqs - omega_plus(n,k)).^2);
        
        % center frequencies
        omega_plus(n+1,k) = (freqs(T/2+1:T)*(abs(u_hat_plus(n+1, T/2+1:T, k)).^2)')/sum(abs(u_hat_plus(n+1,T/2+1:T,k)).^2);
        
    end
    
    % Dual ascent
    lambda_hat(n+1,:) = lambda_hat(n,:) + tau*(sum(u_hat_plus(n+1,:,:),3) - f_hat_plus);
    
    % loop counter
    n = n+1;
    
    % converged yet?
    uDiff = eps;
    for i=1:K
        uDiff = uDiff + 1/T*(u_hat_plus(n,:,i)-u_hat_plus(n-1,:,i))*conj((u_hat_plus(n,:,i)-u_hat_plus(n-1,:,i)))';
    end
    uDiff = abs(uDiff);
     
end
fprintf('%f iterations \n',n);
toc;
time = toc;

%------ Postprocessing and cleanup


% discard empty space if converged early
N = min(N,n);
omega = omega_plus(1:N,:);

% Signal reconstruction
u_hat = zeros(T, K);
u_hat((T/2+1):T,:) = squeeze(u_hat_plus(N,(T/2+1):T,:));
u_hat((T/2+1):-1:2,:) = squeeze(conj(u_hat_plus(N,(T/2+1):T,:)));
u_hat(1,:) = conj(u_hat(end,:));

u = zeros(K,length(t));

for k = 1:K
    u(k,:)=real(ifft(ifftshift(u_hat(:,k))));
end

% remove mirror part
u = u(:,T/4+1:3*T/4);

% recompute spectrum
clear u_hat;
for k = 1:K
    u_hat(:,k)=fftshift(fft(u(k,:)))';
end

end


