% Defining Toeplitz matrices 

% syms s a1 b1;  % Symbolic variable 
% s = sym('s'); 

clear all, close all

f = 0.5;          % Frequency (Hz) Deborah number
w = 2*pi*f;     % Angular frequency
T = 1/f;        % Period


n = 16; % Number of harmonics (-n,...0,...,n)   


nstep = 25;
omega_neg = linspace(-(w-1e-3)/2, -w/(2*nstep), nstep+1); %
omega_pos = linspace(w/(2*nstep), w/2, nstep+1);

omega = [omega_neg 0 omega_pos];




N = 2*n+1; % Matrix dimension 

A = zeros (N,N);
B = zeros (N,N);
I = eye(N);

M = zeros (N,N);


HTM_evals = zeros (N, length (omega));
HTM_evec = cell (length(omega), 1);


for k = 1: N
    M (k,k) = 1i*w*(k-n-1);
end
    

for k = 1: N
    for j = 1: N
        
        if k - j == 0
            A (k,j) = 3/4;
            B (k,j) = 1/8;
        else
        
        idx = k-j;
        c = (exp(-1i*w*idx/2)-1)*(1+0.5*exp(-1i*w*idx/2));
        d = exp(-1i*w*idx/2)*(exp(-1i*w*idx/2)-1);
        
        A(k,j) =  -(1/(1i*w*idx))*c;
        B(k,j) = -(1/(4i*w*idx))*d;
    
        end 
    end
end 

for iter =1: length(omega)

s = omega(iter);

K = (1i*s*I -(A-M));

% Harmonic Transfer Matrix
H = K\I + B;

[evecs, evals] = eig(H);
lambda = diag(evals);

HTM_evals (:, iter) = lambda;
HTM_evec {iter,1} = evecs; 
end
     



    %% Ordering eigenvalues based on proximity

% ordered_evals= evals_ordering (HTM_evals, HTM_evecs); 

k = find(omega == 0);


HTM_evals = sortrows(HTM_evals, k, 'descend' );

ordered_evals= evals_ordering (HTM_evals); 
% ordered_evals= evals_and_evec_ordering (HTM_evals, HTM_evec); 
evals_real = real(ordered_evals);
evals_imag = imag(ordered_evals);


[nrows, ncols] = size(ordered_evals);

idx = (ncols-1)/2+1;

figure()
%Ordered Eigenspectrum
for n =1:nrows
 hold on
    x = evals_real(n,:);
    y = evals_imag(n,:);
    plot(x(:,1:idx-1), y (:, 1:idx-1), '.')
    plot(x (:,idx), y (:, idx), 'o')
    plot(x (:,idx+1:end), y (:, idx+1:end), '.')
    plot(0,0,'+', MarkerSize=10, MarkerFaceColor='r')



end
%  title('Eigenvalue Spectrum - Ordered  nfm = '+string(nfm) +', De = '+string(De)+ ...
% ', \gamma = '+string(gamma0)+ ', Z = '+string(Z))
xlabel('$\mathrm{Re}(\lambda)$', 'Interpreter', 'latex')
ylabel('$\mathrm{Im}(\lambda)$', 'Interpreter', 'latex')
% box on 
% legend ('$\mathrm{\gamma}$ = '+string(gamma0), 'Interpreter', 'latex')
% legend box off