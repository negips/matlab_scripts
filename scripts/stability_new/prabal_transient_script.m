%% Transient Growth-Re

re_range = [200 500 1000 2000 4000];
N_range  = [512 512 1024 2048 4096];

l1 = length(re_range);

for ii = 1:l1
    Re = re_range(ii)
    N = N_range(ii)
    
    run TransientGrowth
    
    save(['trans_growth_re' num2str(Re)]);
    
    clearvars -except re_range ii l1 N_range
    
end