%% Function to generate Random Rate Matrix
% ns = number of states
% K(i,j) rate constant for the i --> j process
% rand is a subroutine that generates a uniformly
% distributed random number between (0,1)
% (a,b) are the states you want to add stochastic bottleneck to

function f = randMat(out,ns,a,b)
P = pwd;

%Output Folder Path
OutFolder = sprintf('%s\\%s', P,out);%
foldCheck(OutFolder);

for i=1:ns
    for j = 1:ns
        if i ~= j 
            K(i,j)=10*rand;
        end
    end
end

% Add a stochastic "bottle neck" between states a and b, by setting smaller rates here
K(a,b)=rand;
K(b,a)=rand;
for i=1:ns
    K(i,i)=-sum(K(:,i));
end

display(K)
f = K;

fout1 = fullfile(OutFolder, sprintf('KMRND%d.txt',ns));
dlmwrite(fout1, K, 'delimiter', '\t');

end
