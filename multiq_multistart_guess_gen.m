function X0s = multiq_multistart_guess_gen(Q,M,n_start_peaks,L,w,lb,ub)
% Generate initial guesses. Generates n_start_peaks evenly spaced peaks and all
% pairs of those peaks.

% Q=size(g2s,1);
% M=length(w)/L;

X = 2*Q+M*L; %length of solution vector

%some peak positions
start_peaks = fix(linspace(1,M,n_start_peaks+2));
start_peaks = start_peaks(2:end-1);
half_width = fix(M/(n_start_peaks+1)/2);
%fix(linspace(4,M-3,n_start_peaks));%4:ceil(M/5):M-3;
%peak pairs
%start_peak_pairs = nchoosek(start_peaks,2);

% n_start_peaks = length(start_peaks);
%n_start_peak_pairs = size(start_peak_pairs,1);
n_start = (1+n_start_peaks)^L;%+n_start_peak_pairs)^L;

one_component_grid = zeros(M,1+n_start_peaks);%+n_start_peak_pairs);

X0s = zeros(X,n_start);
%bl_guesses = zeros(Q, n_start); %set baselines to zero
%bl_guesses = repmat( g2s(:,end), 1, n_start)-1; %set baselines based on
%last data point

%X0s(1:Q,:) = bl_guesses; %baseline guesses
%X0s(Q+1:2*Q,:) = g2s(1,1) - bl_guesses-1; %contrast guesses, use lowest q as reference (most likely to be close to correct)

% contrast and guess centered in bounds. It's better that X0 does not
% depend on the data.
X0s(1:2*Q,:) = repmat((lb(1:2*Q) + ub(1:2*Q))/2,1,n_start);

% impulses
x = (-half_width:half_width);
for i = 1:n_start_peaks
    one_component_grid(start_peaks(i)-half_width:start_peaks(i)+half_width,1+i) = exp(half_width^2./(x.^2-(half_width+1)^2));%0.5;
end
% for i = 1:n_start_peak_pairs
%          one_component_grid(1+start_peak_pairs(i,1)-3:1+start_peak_pairs(i,1)+3,1+n_start_peaks+i) = 0.5;
%          one_component_grid(1+start_peak_pairs(i,2)-3:1+start_peak_pairs(i,2)+3,1+n_start_peaks+i) = 0.5;
% end

% argument list for combvec as cell array
C=cell(L,1);
for l = 1:L
    C{l}=one_component_grid;
end

%combinations per component
X0_combination_grid = combvec(C{:});

X0s(1+2*Q:end,:) = X0_combination_grid;
% for l=1:L
%     for i = 1:n_start_peaks
%         X0s(1+2*Q+M*(l-1)+start_peaks(i)-3:1+2*Q+M*(l-1)+start_peaks(i)+3,1+i) = 0.5;
%     end
%     for i = 1:n_start_peak_pairs
%          X0s(1+2*Q+M*(l-1)+start_peak_pairs(i,1)-3:1+2*Q+M*(l-1)+start_peak_pairs(i,1)+3,1+n_start_peaks+i) = 0.5;
%          X0s(1+2*Q+M*(l-1)+start_peak_pairs(i,2)-3:1+2*Q+M*(l-1)+start_peak_pairs(i,2)+3,1+n_start_peaks+i) = 0.5;
%     end
% end

%normalize guesses
norms=(w*X0s(1+2*Q:end,2:end));
% norms=sum(X0s(1+2*Q:end,2:end));
X0s(1+2*Q:end,2:end)=X0s(1+2*Q:end,2:end)./norms;

%don't return the guess with all zeros.
X0s = X0s(:,2:end);

end
