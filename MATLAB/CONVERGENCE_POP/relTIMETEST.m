%% FUNCTION TO FIND RELAXATION TIME FROM A RATE OR PROB MATRIX
%From the right eigenvector and eigenvalues find the equilibrium
%populations & relaxation times

%in = location of input files (str)
%out = location of output files (str)
%temp = Temperature of simulation (int)
%pep = name of peptide (str)
%data = input Prob/Rate Matrix (str)
%lag = file of lag values needed to convert lag value to time (str)
%tstep = time step of simulation in seconds (int)
%OUTPUTS ARE:   eigenvalues per lag time
%               right & left evectors per lag time
%               relaxation times per lag time
%               vector of slowest rel times
in = 'functionTest\outGILL5';
out ='functionTest\outGILL6';
ns = 2; %number of states
temp = 5; %sim temperature
pep = 'GILL'; %name of peptide

data = sprintf('%d%s_KM*.txt',temp,pep);

totalT = length(data1); %total number of time frames

P = pwd;

% Specify the folder where the files live.
myFolder = sprintf('%s\\%s\\', P, in);
foldCheck(myFolder); %check folder exists

% Specify Output folder
OutFolder =  sprintf('%s\\%s\\', P, out);
foldCheck(OutFolder);

filePattern = fullfile(myFolder, data); % Change to whatever pattern you need.
theFiles = dir(filePattern);
len = length(theFiles); %number of files

lagfile = fullfile(myFolder, lag); %read in lag time step file
lagdata = load(lagfile, '-ascii');

slowRelT = zeros(len , 1); %empty mat to store slowest relT in
%get eigenvalues & populations from the rate matrix


for k = 1 : len
    baseFileName = theFiles(k).name;
    fullFileName = fullfile(myFolder, baseFileName);
    fprintf(1, 'Now reading %s\n', fullFileName);
    % Now do whatever you want with this file name,
    
    M=load(fullFileName,'-ascii') ;
    
    %check to catch any nan or inf values
    for i = 1: ns
        for j = 1:ns
            if isnan(M(i,j)) == 1 || isinf(M(i,j)) ==1
                M(i,j) = 0;
                fileID = fopen(sprintf('%s\\%d%s_errorLog%d.txt',OutFolder, temp,pep,k), 'w');
                formatSpec = '%d_TM(%d,%d) = NaN/Inf\n';
                fprintf(fileID,formatSpec,k,i,j);
            end
        end
    end
    
    
    [V, D] = eig(M); %V = right evectors, D = eigenvalues
    
    %the Prob Matrix is a left stochastic matrix
    %square real matrix with columns summing to 1
    %it is not symmetric and V and W will not nec be the same
    
    %convert it to time by multiplying by the t step.
    
    %sort eigenvalues and vectors in descending order s.t. eval = 1 is first
    [D,idx] = sort(diag(D), 'descend');
    V = V(:, idx);
    
   
    rel_time = -1*(lagdata(k)./D); % work out relaxation time in frames
    slowRelT(k) = rel_time(2);
    
    % 1/lamda = relaxation time, second largest evalue gives slowest rel time
    %one evalue will be 1, rest will be <1 for PM
    % one eval will be 0, rest <1 for KM
    
    fout1 = fullfile(OutFolder, sprintf('%d%s_EVAL%06d.txt',temp,pep,k));
    fout2 = fullfile(OutFolder, sprintf('%d%s_REVEC%06d.txt',temp,pep,k));
    fout4 = fullfile(OutFolder, sprintf('%d%s_RELT%06d.txt',temp,pep,k));
    
    dlmwrite(fout1, V, 'delimiter', '\t');
    dlmwrite(fout2, D, 'delimiter', '\t');
    dlmwrite(fout4, rel_time, 'delimiter', '\t');

end
    %write out slowest relaxation time vector
    fout5 = fullfile(OutFolder,sprintf('%d%s_RT.txt',temp,pep));
    dlmwrite(fout5, slowRelT, 'delimiter', '\t');
    f = [lagdata' , slowRelT];
    
    fout6 = fullfile(OutFolder,sprintf('%d%s_LTRT.txt',temp,pep));
    dlmwrite(fout6, f, 'delimiter','\t');
    


%