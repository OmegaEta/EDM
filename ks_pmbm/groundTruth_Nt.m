
N_t = zeros(100,1);
for t=1:100
    
    N_t(t) = size(groundTruth{t}.x,2);
end

