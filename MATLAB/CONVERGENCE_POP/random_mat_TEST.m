clear variables
ns = 3;

for i=1:ns
    for j = 1:ns
        if i ~= j 
            K(i,j)=10*rand;
        end
    end
end

for i=1:ns
    K(i,i)=-sum(K(:,i));
end

disp(K)
