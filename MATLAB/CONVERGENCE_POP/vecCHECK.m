
clear variables
iter = 0;
lenMat = round(9/2);
c = zeros(lenMat,1);
for lag = 1:2:9

        iter = iter+1;
        c(iter) =lag; %save out lag values
        
end
