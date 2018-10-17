
%get eigenvalues & populations from the symmetrised count matrix
f_in='dcsa_510_CM_SYM_1.txt';
CM_data=load(f_in,'-ascii') ; 

%variables:
temp = 510; %temperature of simulation
pep = 'dcsa'; %csa or dcsa peptide
lag = 1; %lag time
tstep = 500e-12; %time step

CM_eigs = eig(CM_data);

time = lag*tstep; %total time

rel_t = CM_eigs.*time; % work out relaxation time 

Pop = sum(CM_data); % gives population of each state

file_name1 = sprintf('%d_%s_CM_eigs_%d.txt',temp,pep,lag);
file_name2 = sprintf('%d_%s_CM_Pop_%d.txt',temp,pep,lag);
file_name3 = sprintf('%d_%s_CM_rel_t_%d.txt',temp,pep,lag);

dlmwrite(file_name1, CM_eigs, 'delimiter', '\t'); 

dlmwrite(file_name2, Pop, 'delimiter', '\t'); 

dlmwrite(file_name3, rel_t, 'delimiter', '\t'); 
