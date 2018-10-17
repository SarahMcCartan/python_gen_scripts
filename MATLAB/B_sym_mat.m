close all


%read in count matrix created after TBA applied and symmetrize the matrix

f_in='310_Acsa_count_matrix_TBA_1.txt';
CM=load(f_in,'-ascii') ; 

%variables per simulation
ns = 3; %number of states
temp = 310; %temperature of simulation
pep = 'Acsa'; %csa or dcsa peptide
lag = 1;


POP = sum(CM);

CM_sym=zeros(ns,ns);

%calc averages of C(ij) and (Cji)


for i=1:ns
    for j=1:ns
        if i ~= j
            
            CM_sym(i,j) = (CM(i,j) + CM(j,i))/2 ; 
            
        end
        
    end
    
end



%Diagonal of new matrix to be same as old matrix

for x = 1:ns
    CM_sym(x,x) = CM(x,x);
    
end

POP2 = sum(CM_sym); %sum of each col

%adjust the diagonal count to keep total pop as unsymm matrix

for y = 1:ns
    if POP(y) ~= POP2(y)
        CM_sym(y,y) = CM_sym(y,y) + (POP(y)-POP2(y));
        y;
    end
    
end

POP3 = sum(CM_sym); %check the new sum is as per POP

if POP == POP3
    disp('populations consistent')
    
end

file_name1 = sprintf('%d_%s_CM_SYM_%d.txt',temp,pep,lag);

dlmwrite(file_name1, CM_sym, 'delimiter', '\t'); 


