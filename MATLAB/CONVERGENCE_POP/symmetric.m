
%% FUNCTION TO MAKE A SYMMETRIC COUNT MATRIX
% THIS IS TO SATISFY DETAILED BALANCE
%i.e.  Cij = Cji BY (Cij+Cji)/2
%in = location of input files
%out = location of output files
%temp = Temperature of simulation
%pep = name of peptide
%ns = number of states
%data = input unsymmetrised Count Matrix
%OUTPUTS ARE:   sym = Symmetrised Matrix of No. of Transitions per lag value

function s = symmetric(in,out,data,ns,temp,pep)

P = pwd;

% Specify the folder where the files live.
myFolder = sprintf('%s\\%s\\', P,in);
% Check to make sure that folder  exists
foldCheck(myFolder);

OutFolder = sprintf('%s\\%s', P, out);
% Check to make sure that folder exists
foldCheck(OutFolder);

%Get a list of all files in the folder with the desired file name pattern.
filePattern = fullfile(myFolder, data); % Change to whatever pattern you need.
theFiles = dir(filePattern);


for k = 1 : length(theFiles)
    baseFileName = theFiles(k).name;
    fullFileName = fullfile(myFolder, baseFileName);
    fprintf(1, 'Now reading %s\n', fullFileName);
    % Now do whatever you want with this file name,
    % such as reading it in as an image array with imread()
    
    CM=load(fullFileName,'-ascii') ;
    
    %POP = sum(CM); %population is the sum of the cols
    CM_sym=zeros(ns,ns); %create empyty matrix
    
    %calc averages of C(ij) and (Cji)
    
    for i=1:ns
        for j=1:ns
            if i ~= j
                CM_sym(i,j) = (CM(i,j) + CM(j,i))/2 ;
                
            end
            
        end
        
    end
    
    %Diagonal of new matrix to be same as old matrix
    
    for x = 1:ns
        CM_sym(x,x) = CM(x,x);
    end
    %
    %     POP2 = sum(CM_sym); %sum of each col
    %
    %     %adjust the diagonal count to keep total pop as unsymm matrix
    %
    %     for y = 1:ns
    %       if POP(y) ~= POP2(y)
    %             CM_sym(y,y) = CM_sym(y,y) + (POP(y)-POP2(y));
    %        end
    %     end
    %
    %     POP3 = sum(CM_sym); %check the new sum is as per original POP
    %
    %     if POP == POP3
    %         disp('Pop OK')
    %     else
    %         disp('Pop Error')
    %     end
    
    s = CM_sym;
    
    fout1 = fullfile(OutFolder, sprintf('%d%s_CMSYM%06d.txt',temp,pep,k));
    dlmwrite(fout1, CM_sym, 'delimiter', '\t');
    
    
end


end

