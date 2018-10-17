 
in = 'functionTest\outGILL7\CMS';
in2 = 'functionTest\outGILL7\CM';
out ='functionTest\outGILL7\KM3';

ns = 2; %number of states
temp = 7; %sim temperature
pep = 'GILL'; %name of peptide
lag = sprintf('%d%s_lagVals.txt',temp,pep);


data = sprintf('%d%s_CMS*.txt',temp,pep);

P = pwd; %find current dir

% Specify the folder where the files live.
myFolder =  sprintf('%s\\%s\\', P,in);
foldCheck(myFolder); %check folder exists

% Specify the folder where the files live.
myFolder2 =  sprintf('%s\\%s\\', P,in2);
foldCheck(myFolder2); %check folder exists

%specify output folder
OutFolder =  sprintf('%s\\%s', P,out);
foldCheck(OutFolder);

%Get a list of all files in the folder with the desired file name pattern.
filePattern = fullfile(myFolder, data); % Change to whatever pattern you need.
theFiles = dir(filePattern);

lagfile = fullfile(myFolder2, lag); %read in lag time step file
lagdata = load(lagfile, '-ascii');

for k = 1 : length(theFiles)
    baseFileName = theFiles(k).name;
    fullFileName = fullfile(myFolder, baseFileName);
    fprintf(1, 'Now reading %s\n', fullFileName);
    % Now do whatever you want with this file name,
    
    %PM=load(fullFileName,'-ascii') ; %load in Prob Mat
    CMS = load(fullFileName,'-ascii');
    tau = lagdata(k); %relevant lag time
    
    %KM = (-1/tau).*(log(1-PM));
    KM =  (-1/tau).*(log(CMS));
    for x = 1:ns
        for y = 1:ns
            if x == y
                KM(x,x) = 0; %remove diagonal values
            end
        end
    end
%     
    KSUM2 = sum(KM); %get totals for each col in KM
    
    for z = 1:ns
        KM(z,z) = 0 - KSUM2(z); % all columns must sum to 0
    end
%     
%     for w = 1:ns
%         
%         if sum(K(:,w)) == 0
%             disp(w)
%             disp('KM cols sum to 0')
%         else
%             disp(w)
%             disp('KM Error')
%         end
%         
%     end

    fout2 = fullfile(OutFolder,sprintf('%d%s_KM%05d.txt',temp,pep,k));
    dlmwrite(fout2, KM, 'delimiter', '\t');
    
    
end


fprintf('Step 4: Find Relaxation Time\n')

data4 = sprintf('%d%s_KM*.txt',temp,pep);

relTime(out,out,data4,lag,ns,temp,pep);

fprintf('Step 5: Plot Slowest Rel Time vs Lag Time to check Convergence\n')

data5 = sprintf('%d%s_LTRT.txt', temp,pep);

relConv(out,out,data5,tstep,temp,pep);

fclose('all'); %close any open error log files
    