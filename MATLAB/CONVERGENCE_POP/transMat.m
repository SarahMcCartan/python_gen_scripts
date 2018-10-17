%% FUNCTION TO BUILD A TRANSITION PROB MATRIX FROM A SYMMETRIC COUNT MATRIX
%in = location of input files
%out = location of output files
%temp = Temperature of simulation
%pep = name of peptide
%ns = number of states
%data = input Symmetrised Count Matrix
%OUTPUTS ARE:   TM(TAU) = TRANSITION PROB matrix per lag val (TAU)


function TM = transMat(in,out,data,ns,temp,pep)

P = pwd; %find current dir

% Specify the folder where the files live.
myFolder =  sprintf('%s\\%s\\', P,in);
foldCheck(myFolder); %check folder exists

%specify output folder
OutFolder =  sprintf('%s\\%s', P,out);
foldCheck(OutFolder);

%Get a list of all files in the folder with the desired file name pattern.
filePattern = fullfile(myFolder, data); % input files = CM sym
theFiles = dir(filePattern);


for k = 1 : length(theFiles)
    baseFileName = theFiles(k).name;
    fullFileName = fullfile(myFolder, baseFileName);
    fprintf(1, 'Now reading %s\n', fullFileName);
    % Now do whatever you want with this file name,
    
    CMS=load(fullFileName,'-ascii') ; %load in Count matrix symmetric
    TM = zeros(ns,ns);
    
    for i =1:ns
        for j =1:ns
            TM(i,j) = CMS(i,j)/(sum(CMS(:,j)));
        end
    end
   
    fout1 = fullfile(OutFolder,sprintf('%d%s_TM%06d.txt',temp,pep,k));
    dlmwrite(fout1, TM, 'delimiter', '\t');
    
    
end

end




