clear variables

in = '310_CSA_CG\KinAnalysisTESTKM';
out ='310_CSA_CG\KinAnalysisTESTKM';
ns = 2; %number of states
temp = 310; %sim temperature
pep = 'CSA';
data = sprintf('%d%s_KA*.txt',temp,pep);
lag = sprintf('%d%s_lagVals.txt',temp,pep);
tstep = 500e-12;

P = pwd; %find current dir

% Specify the folder where the files live.
myFolder = sprintf('%s\\%s', P, in);
foldCheck(myFolder); %check folder exists

% Specify output folder
OutFolder = sprintf('%s\\%s', P, out);
foldCheck(OutFolder);



lagFile = fullfile(myFolder, lag); %read in lag val file
lagData = load(lagFile, '-ascii');

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
    
    %infNan(OutFolder,ns,M,k,temp,pep);
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
    
    
    [V, D, W] = eig(M); %V = right evectors, D = eigenvalues, W = left evectors
    
    %the Prob Matrix is a left stochastic matrix
    %square real matrix with columns summing to 1
    %it is not symmetric and V and W will not nec be the same
    
    time = lagdata(k)*tstep; %load the correct lag time for the TM created
    %convert it to time by multiplying by the t step.
    
    %sort eigenvalues and vectors in descending order s.t. eval = 1 is first
    [D,idx] = sort(diag(D), 'descend');
    V = V(:, idx);
    %W = W(:,idx);
    
    %[eigvec,eigval]=eig(M); % diagonalize K, eigvec stores the eigenvectors, eigval the eigenvalues

    %[dsorted,index]=sort(diag(eigval),'descend'); % sort the eigenvalues. dsorted stores the eigenvalues, index the corresponding indices
    % sorted eigenvalues:
    %ind=index(1); %index corresponding to 0 eigenvalue
    %eq=eigvec(:,ind)/sum(eigvec(:,ind));

    %[eigvec,eigval]=eig(M'); % diagonalize K, eigvec stores the right eigenvectors
    %slowest_relrate=-dsorted(2);
    %slow_vec=eigvec(:,index(2));

   
    rel_time = time*(1./D); % work out relaxation time
    slowRelT(k) = rel_time(2);
    % 1/lamda = relaxation time, second largest evalue gives slowest rel time
    %one evalue will be 1, rest will be <1 for PM
    % one eval will be 0, rest <1 for KM
    
    %fout1 = fullfile(OutFolder, sprintf('%d%s_EVAL%06d.txt',temp,pep,k));
    %fout2 = fullfile(OutFolder, sprintf('%d%s_REVEC%06d.txt',temp,pep,k));
    %fout3 = fullfile(OutFolder, sprintf('%d%s_LEVEC%06d.txt',temp,pep,k));
    %fout4 = fullfile(OutFolder, sprintf('%d%s_RELT%06d.txt',temp,pep,k));
    
    %dlmwrite(fout1, V, 'delimiter', '\t');
    %dlmwrite(fout2, D, 'delimiter', '\t');
    %dlmwrite(fout3, W, 'delimiter', '\t');
    %dlmwrite(fout4, rel_time, 'delimiter', '\t');

end
    %write out slowest relaxation time vector
    %fout5 = fullfile(OutFolder,sprintf('%d%s_RT.txt',temp,pep));
    %dlmwrite(fout5, slowRelT, 'delimiter', '\t');
    f = [lagdata' , slowRelT];
    %fout6 = fullfile(OutFolder,sprintf('%d%s_LTRT.txt',temp,pep));
    %dlmwrite(fout6, f, 'delimiter','\t');