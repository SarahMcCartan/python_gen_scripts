%% Script for Kinetic Analysis of MD trajectories
close all
clear variables

%%
T = 310;
name = 'CSA';
fprintf('Step 1: Count Matrix by Sliding Window technique\n')

count('functionTest','functionTest\outCSA310_3',310,'CSA',2,'310_CSA_CG_TRAJ.txt',10,10,1500,1,1,0);
%%
fprintf(' Step 2: Symmetrise CM\n')
%from now on input file is same as output file as using files from pervious
%step
symmetric('functionTest\outCSA310_3','functionTest\outCSA310_3',sprintf('%d%s_CM*.txt',T,name),2,310,'CSA');

%%
fprintf('Step 3: Construct Rate Matrix \n')


rate('functionTest\outCSA310_3','functionTest\outCSA310_3',sprintf('%d%s_CMSYM*.txt',T,name),'310_CSA_CG_TRAJ.txt',2,310,'CSA');

%%
fprintf('Step 4: Find Relaxation Time\n')

relTime('functionTest\outCSA310_3','functionTest\outCSA310_3',sprintf('%d%s_KM0*.txt',T,name),'310CSA_lagVals.txt',2,310,'CSA');

%%
fprintf('Step 5: Plot Slowest Rel Time vs Lag Time to check Convergence\n')

relConv('functionTest\outCSA310_3','functionTest\outCSA310_3',sprintf('%d%s_LTRT.txt',T,name),1,310,'CSA');

fclose('all'); %close any open error log files
%%
fprintf('Kinetic Analysis Finished\n')





