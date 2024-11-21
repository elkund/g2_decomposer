function X0s = multiq_multistart_guess_gen(g2s,L,w)
% Generate initial guesses. Generates some evenly spaced peaks and all
% pairs of those peaks.

Q=size(g2s,1);
M=length(w)/L;

X = 2*Q+M*L; %length of solution vector

%some peak positions
start_peaks = 4:ceil(M/6):M-3;
%peak pairs
start_peak_pairs = nchoosek(start_peaks,2);

n_start_peaks = length(start_peaks);
n_start_peak_pairs = size(start_peak_pairs,1);
n_start = 1+n_start_peaks+n_start_peak_pairs;

X0s = zeros(X,n_start);
bl_guesses = zeros(Q, n_start); %set baselines to zero
%bl_guesses = repmat( g2s(:,end), 1, n_start)-1; %set baselines based on
%last data point

X0s(1:Q,:) = bl_guesses; %baseline guesses
X0s(Q+1:2*Q,:) = g2s(1,1) - bl_guesses-1; %contrast guesses, use lowest q as reference (most likely to be close to correct)

for l=1:L
    for i = 1:n_start_peaks
        X0s(1+2*Q+M*(l-1)+start_peaks(i)-3:1+2*Q+M*(l-1)+start_peaks(i)+3,1+i) = 0.5;
    end
    for i = 1:n_start_peak_pairs
         X0s(1+2*Q+M*(l-1)+start_peak_pairs(i,1)-3:1+2*Q+M*(l-1)+start_peak_pairs(i,1)+3,1+n_start_peaks+i) = 0.5;
         X0s(1+2*Q+M*(l-1)+start_peak_pairs(i,2)-3:1+2*Q+M*(l-1)+start_peak_pairs(i,2)+3,1+n_start_peaks+i) = 0.5;
    end
end

%normalize guesses
norms=(w*X0s(1+2*Q:end,2:end));
X0s(1+2*Q:end,2:end)=X0s(1+2*Q:end,2:end)./norms;

end