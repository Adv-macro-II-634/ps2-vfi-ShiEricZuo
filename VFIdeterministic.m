close all
clear all
%%%% Set up parameters
alpha = 0.35;
beta = 0.99;
delta = 0.025;
sigma = 2;
a_high = 1.1;
a_low = 0.678; %questionable
prob = [0.977 0.023; 0.074 0.926];

%%%% Set up discretized state space
k_min = 0;
k_max = 45;
num_k = 1000; % number of points in the grid for k

k = linspace(k_min, k_max, num_k);

k_mat = repmat(k', [1 num_k]); % this will be useful in a bit

%%%% Set up consumption and return function
% 1st dim(rows): k today, 2nd dim (cols): k' chosen for tomorrow
cons_low = a_low*k_mat' .^ alpha + (1 - delta) * k_mat' - k_mat; 
cons_high = a_high*k_mat' .^ alpha + (1 - delta) * k_mat' - k_mat; 

cons_low(find(cons_low<=0)) = NaN;
cons_high(find(cons_high<=0)) = NaN;

ret_low = cons_low .^ (1 - sigma) / (1 - sigma);
ret_high = cons_high .^ (1 - sigma) / (1 - sigma);% return function
% negative consumption is not possible -> make it irrelevant by assigning
% it very large negative utility
ret_low(find(isnan(ret_low))) = -inf;
ret_high(find(isnan(ret_high))) = -inf;

%%%% Iteration
dis = 1; tol = 1e-07; % tolerance for stopping 
v_guess = zeros(num_k,2);
pol_index = zeros(num_k,2);
while dis > tol
    % compute the utility value for all possible combinations of k and k':
     [vfn_low,pol_index_low]=max(ret_low + beta*repmat(v_guess*prob(2,:)',1,num_k));
  [vfn_high,pol_index_high]=max(ret_high + beta*repmat(v_guess*prob(1,:)',1,num_k));
    % find the optimal k' for every k:
 
  
  pol_index=[pol_index_low' pol_index_high'];
 vfn=[vfn_low' vfn_high'];
  
    
    % what is the distance between current guess and value function
    dis = max(abs(vfn - v_guess));
    
    % if distance is larger than tolerance, update current guess and
    % continue, otherwise exit the loop
    v_guess = vfn;
end

%g = k(pol_indx); % policy function

plot(k,vfn_low,'-',k,vfn_high,':')
figure
plot(k,pol_index_low,'-',k,pol_index_high,':')


