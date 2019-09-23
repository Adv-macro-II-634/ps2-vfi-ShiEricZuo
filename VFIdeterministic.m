close all
clear all
%%%% Set up parameters
alpha = 0.35;
beta = 0.99;
delta = 0.025;
sigma = 2;
a_high = 1.1;
a_low = 0.678;
prob = [0.977 0.023; 0.074 0.926];

%%%% Set up discretized state space
k_min = 0;
k_max = 45;
num_k = 1000; % number of points in the grid for k

k = linspace(k_min, k_max, num_k);

k_mat = repmat(k', [1 num_k]); % this will be useful in a bit

%%%% Set up consumption and return function
% 1st dim(rows): k today, 2nd dim (cols): k' chosen for tomorrow
cons_low = a_low*(k_mat .^ alpha) + (1 - delta) * k_mat - k_mat'; 
cons_high = a_high*(k_mat .^ alpha) + (1 - delta) * k_mat - k_mat'; 

%cons_low(find(cons_low<=0)) = NaN;
%cons_high(find(cons_high<=0)) = NaN;

ret_low = cons_low .^ (1 - sigma) / (1 - sigma);
ret_high = cons_high .^ (1 - sigma) / (1 - sigma);% return function
% negative consumption is not possible -> make it irrelevant by assigning
% it very large negative utility
ret_low(cons_low<0) = -inf;
ret_high(cons_high<0) = -inf;

%%%% Iteration
dis = 1; tol = 1e-06; % tolerance for stopping 
v_guess = ones(2,num_k);
pol_index = ones(2,num_k);
while dis > tol
    % compute the utility value for all possible combinations of k and k':
    vfn_low_mat = ret_low + beta*repmat(prob(2,:)*v_guess,[num_k,1]);
    vfn_high_mat = ret_high + beta*repmat(prob(1,:)*v_guess,[num_k,1]);
    
    [vfn_low,pol_index_low]=max(vfn_low_mat,[],2);
  [vfn_high,pol_index_high]=max(vfn_high_mat,[],2);
    % find the optimal k' for every k:
 
  
  pol_index=[pol_index_low'; pol_index_high'];
 vfn=[vfn_high';vfn_low'];
  
    
    % what is the distance between current guess and value function
    dis_mat=abs(vfn - v_guess);
    dis = max(dis_mat(:));
    
    % if distance is larger than tolerance, update current guess and
    % continue, otherwise exit the loop
    v_guess = vfn;
end

g_h = k(pol_index_high); % policy function
g_l = k(pol_index_low);
plot(k,vfn_low,'-',k,vfn_high,':')
figure
plot(k,g_l,'-',k,g_h,':')

saving_h = g_h-(1-delta).*k;
saving_l = g_l-(1-delta).*k;

figure
plot(k,saving_l,'-',k,saving_h,':')


%simulation
T_sim=1000;
a_high_sim=1.0001;
a_low_sim=(1-a_high_sim*0.7629)/0.2371;
a_state=[a_high_sim;a_low_sim];
%draw random number for A
rng(1);
rand_nums=rand(T_sim,1);

%turn random numbers into value for A using 
A_sim = zeros(T_sim,1);
A_sim(1)=1;
%start with arbitrary capital stock, then follow the policy funtion
%according to simulated state in the current period
k_sim=zeros(T_sim,1);
k_sim(1)=5;

%Space to store the simulated values of tech and capital
A=zeros(T_sim,1);
K=zeros(T_sim,1);

%Process of simulation

for t=1:T_sim
    if A_sim(t) ==1
        if  rand_nums(t)<0.977 
               A_sim(t+1) = 1; 
        else 
               A_sim(t+1) = 2; 
        end
    elseif rand_nums(t)<0.926 
               A_sim(t+1) = 2; 
             else 
                 A_sim(t+1) =1;
    end
  k_sim(t+1)=pol_index(A_sim(t),k_sim(t));            
end


%Update the simulated values into the space
A=a_state(A_sim);
K=k(k_sim);

%Calculate output 
y=A'.*(K.^alpha);

%remove the first few, unusually small y
y(1:20)=[];

%Calculate s.d.
sdy=std(y);



