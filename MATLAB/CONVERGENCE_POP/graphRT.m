

%% FUNCTION TO PLOT REL T VS LAG T TO CHECK FOR CONVERGENCE OF REL T
%in = location of input files (str)
%out = location of output files (str)
%temp = Temperature of simulation (int)
%pep = name of peptide (str)
%data = lag time vs rel time matrix
%tstep = time step of simulation in seconds (int)
%OUTPUTS ARE:   lag vs rel time text file
%               graph
in = 'functionTest\outGILL3';
out ='functionTest\outGILL3';
ns = 2; %number of states
temp = 1; %sim temperature
pep = 'GILL'; %name of peptide
tstep = 1;
data = sprintf('%d%s_LTRT.txt', temp,pep);

iter = 2;
P = pwd; %find current dir

% Specify the folder where the files live.
myFolder = sprintf('%s\\%s', P, in);
foldCheck(myFolder); %check folder exists

% Specify output folder
OutFolder = sprintf('%s\\%s', P, out);
foldCheck(OutFolder);

lagRelIn = fullfile(myFolder, data); %read in lag val file
lagRelData = load(lagRelIn, '-ascii');

lagRelTime = tstep.*(lagRelData);
%graph to check convergence of relaxation time for long lag time
figure()

x = lagRelTime((1:500),1);
y = lagRelTime((1:500),2);

plot(x , y)

xlabel('\tau (s)')
ylabel('\mu (s)')

fout2 = fullfile(OutFolder,sprintf('%d%s_RelvsLagT%d.png', temp,pep,iter));
saveas(gcf,fout2 );

fout3 = fullfile(OutFolder,sprintf('%d%s_RelvsLagT%d.fig', temp,pep,iter));
saveas(gcf,fout3 );

