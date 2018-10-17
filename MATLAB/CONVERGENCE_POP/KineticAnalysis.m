%% Script for Kinetic Analysis of MD trajectories
close all
clear variables

% Step 1: Count Matrix by Sliding Window technique

in = 'functionTest';
out ='functionTest\outSarah';
ns = 2; %number of states
temp = 310; %sim temperature
pep = 'SAZZ'; %name of peptide

frame = 1; 
a=1;
b=5;
c=10;
sd =1;
data1 = 'SARAH.txt'; %TBA trajectory

S1 = count(in,out,temp,pep,ns,data1,a,b,c,frame,sd);
disp(S1)

% Step 2: Symmetrise CM
%from now on input file is same as output file as using files from pervious
%step

data2 = sprintf('%d%sCM*.txt',temp,pep);

S2=symmetric(out,out,data2,ns,temp,pep);
disp(S2)

% Step 3: Construct Prob Matrix and Rate Matrix

data3 =  sprintf('%d%sCMSYM*.txt',temp,pep);
lag = sprintf('%d%slagVals.txt',temp,pep);

probmat(out,out,data3,lag,ns,temp,pep);

% Step 4: Find Relaxation Time

data4 = sprintf('%d%sPM*.txt',temp,pep);
tstep = 500e-12; 

relTime(out,out,data4,lag,tstep,ns,temp,pep);

% Step 5: Plot Slowest Rel Time vs Lag Time to check Convergence

data5 = sprintf('%d%sRT.txt', temp,pep);

relConv(out,out,data5,lag,tstep,temp,pep);









