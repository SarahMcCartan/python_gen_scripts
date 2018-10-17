

%% FUNCTION TO PLOT REL T VS LAG T TO CHECK FOR CONVERGENCE OF REL T
%in = location of input files (str)
%out = location of output files (str)
%temp = Temperature of simulation (int)
%pep = name of peptide (str)
%data = lag time vs rel time matrix
%tstep = time step of simulation in seconds (int)
%OUTPUTS ARE:   lag vs rel time text file
%               graph

function rc = relConv(in,out,data,tstep, temp, pep)

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

x = lagRelTime(:,1);
y = lagRelTime(:,2);

plot(x , y)

xlabel('\tau (s)')
ylabel('\mu (s)')

fout2 = fullfile(OutFolder,sprintf('%d%s_RelvsLagT.png', temp,pep));
saveas(gcf,fout2 );

fout3 = fullfile(OutFolder,sprintf('%d%s_RelvsLagT.fig', temp,pep));
saveas(gcf,fout3 );

rc = lagRelTime;


end

