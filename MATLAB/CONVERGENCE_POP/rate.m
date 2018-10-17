%% FUNCTION TO BUILD A RATE MATRIX FROM A SYMMETRIC TRANSITION COUNT MATRIX
%in = location of input files
%out = location of output files
%temp = Temperature of simulation
%pep = name of peptide
%ns = number of states
%trajdata = trajectory txt file
%data = input Symmetrised Count Matrix
%OUTPUTS ARE:   KM = rate matrix per lag val

function km = rate(in,out,data,traj,ns,temp,pep)

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

trajfile = fullfile(myFolder, traj); %read in lag time step file
trajdata = load(trajfile, '-ascii');
tt = length(trajdata);

for k = 1 : length(theFiles)
    baseFileName = theFiles(k).name;
    fullFileName = fullfile(myFolder, baseFileName);
    fprintf(1, 'Now reading %s\n', fullFileName);
    % Now do whatever you want with this file name,
    
    CMS=load(fullFileName,'-ascii') ; %load in Count matrix symmetric
    total = sum(sum(CMS)); %total time in frames
    
    %Find Peq estimate
    Peq = zeros(ns,1);
    for q = 1:ns
        Peq(q) = sum(CMS(:,q))./total;
    end
    
    KM = zeros(ns,ns);
    %KM(i,j) = CMS(i,j)./(Peq(i)*TotalT)
    
    for i =1:ns
        for j =1:ns
            if i ~= j
                KM(i,j) = rdivide(CMS(i,j),(Peq(j).*(tt)));
            end
        end
    end
    
    KSUM = sum(KM); %get totals for each col in PM
    
    for p = 1:ns
        KM(p,p) = 0 - KSUM(p); % all columns should sum to 0
    end
    
    %check sums of each column equal 1:
    
    for r = 1:ns
        
        if sum(KM(:,r)) == 0
            disp(r)
            disp('KM cols sum to 0')
        else
            disp(r)
            disp('KM Error')
        end
        
    end
    
    fout1 = fullfile(OutFolder,sprintf('%d%s_KM%06d.txt',temp,pep,k));
    dlmwrite(fout1, KM, 'delimiter', '\t');
    
    fout2 = fullfile(OutFolder,sprintf('%d%s_PeqEST%06d.txt',temp,pep,k));
    dlmwrite(fout2, Peq, 'delimiter', '\t');
    
   km = KM;
   
end

end




