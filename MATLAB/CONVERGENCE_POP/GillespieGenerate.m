clear variables

%Script to generate a Gillespie Trajectory

in = 'functionTest';
out ='functionTest';
ns = 2; %number of states
temp = 2; %sim temperature
pep = ''; %name of peptide

NN = 100; %no. of transitions
bins =100; %for tDisc

% generate random matrix as a test
%f = randMat('functionTest',ns,1,2);
%data1 = f;

data1 = 'KMRND2.txt'; %rate matrix either real or generated

fprintf 'Generate Gillespie Trajectory\n'

g = gillTraj(in,out,data1,ns,NN,temp,pep);

fprintf 'Discretize Time Step\n'

data2 = sprintf('%d%s_gillTraj_rough.txt',temp,pep); 

td = tDisc(out,out,data2,bins,NN,temp,pep);

fprintf 'Finished Gillespie Generate'


