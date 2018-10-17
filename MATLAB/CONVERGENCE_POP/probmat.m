%% FUNCTION TO BUILD A PROBABILITY MATRIX FROM A SYMMETRIC TRANSITION COUNT MATRIX
%in = location of input files
%out = location of output files
%temp = Temperature of simulation
%pep = name of peptide
%ns = number of states
%data = input Symmetrised Count Matrix
%OUTPUTS ARE:   PM = prob matrix per lag val

function pm = probmat(in,out,data,ns,temp,pep)

P = pwd; %find current dir

% Specify the folder where the files live.
myFolder =  sprintf('%s\\%s\\', P,in);
foldCheck(myFolder); %check folder exists

%specify output folder
OutFolder =  sprintf('%s\\%s', P,out);
foldCheck(OutFolder);

%Get a list of all files in the folder with the desired file name pattern.
filePattern = fullfile(myFolder, data); % Change to whatever pattern you need.
theFiles = dir(filePattern);




for k = 1 : length(theFiles)
    baseFileName = theFiles(k).name;
    fullFileName = fullfile(myFolder, baseFileName);
    fprintf(1, 'Now reading %s\n', fullFileName);
    % Now do whatever you want with this file name,
    
    CMS=load(fullFileName,'-ascii') ; %load in Count matrix symmetric
    
    
    %% PROB MATRIX CONSTRUCT
    PM = zeros(ns,ns);
    
    for i =1:ns
        for j =1:ns
            if i ~= j
                PM(i,j) = rdivide(CMS(i,j),sum(CMS(:,j)));
            end
        end
    end
    
    PSUM = sum(PM); %get totals for each col in PM
    
    for p = 1:ns
        PM(p,p) = 1 - PSUM(p); % all columns should sum to 1
    end
    
    %check sums of each column equal 1:
    
    for q = 1:ns
        
        if sum(PM(:,q)) == 1
            disp(q)
            disp('PM cols sum to 1')
        else
            disp(q)
            disp('PM Error')
        end
        
    end
    
    fout1 = fullfile(OutFolder,sprintf('%d%s_PM%06d.txt',temp,pep,k));
    dlmwrite(fout1, PM, 'delimiter', '\t');
    
   pm = PM;
end

end

