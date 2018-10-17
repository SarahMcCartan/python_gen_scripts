close all

%MASTER MATLAB CODE 2 
%RATE MATRIX GENERATION

f_in='310_Acsa_CM_SYM_1.txt';
f_in2 ='310_Acsa_Life_T_Avg_1.txt';
CMS=load(f_in,'-ascii') ; 
Life_T = load(f_in2, '-ascii');

%variables
temp = 510; %temperature of simulation
pep = 'dcsa'; %csa or dcsa peptide
lag = 1; %lag time
%LB method uses lag time, have defined the lagtime here, set to 1 for now.
%define lagtime in n steps which is n*500ps (time between each step)
R = 0.15;
ns = 3;

N=sum(CMS); %total number of transitions from l to i where l goes from l to
%3

%gives transition matrix where PM(i,j)= CM(i,j)/sumC(i,l)
%here sum(C(i,l))==N(i)


PM=zeros(ns,ns);

for i = 1:ns
   for j = 1:ns
       if i == j

        PM(i,j)=CMS(i,j)./(N(i));
       else
        PM(i,j) = CMS(i,j)./N(i);
        
       end
   end
end


fname13 =  sprintf('%d_%s_PM_%d.txt',temp,pep,lag);

dlmwrite(fname13, PM, 'delimiter', '\t'); %write out to txt file


%create rate matrix K
K=zeros(ns,ns);

for k = 1:ns
    for l = 1:ns
        if k ==l
            K(k,l) = -1./(Life_T(k));
        else
            K(k,l) = CMS(k,l)./N(k);
        end
    end
end


%K(i,j) = BP(i,j)/T(i)
%branching probability of going from state i to j
%is N(i,j)sym/(sum of N(l,i)) where N is number of transitions
%the below gives output of a mathematically correct rate matrix or
%infinetismal generator, however I have not divided by the lifetime..need
%to investigate more


fname14 =  sprintf('%d_%s_Rate_Matrix_%d.txt',temp,pep,lag);
dlmwrite(fname14, K, 'delimiter', '\t'); %write out to txt file

%%
%find stationary distribution of K which is the equilibrium pop
%singular value decompositing to get peq that satisifies K*peq=0
%this is not working, need to get non trivial answer here, matlab returning
% all zeros. Can work out by hand but not ideal

%[U S V] = svd(inv(K));
%peq=V(:,end);

