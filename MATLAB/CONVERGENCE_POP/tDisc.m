%% FUNCTION TO DISCRETIZE TIME STEP OF A TRAJ
%Code to create constant time step trajectory from the Gillespie traj
%which contains the length of time each state lived for before
%transitioning to new state
%%in = location of input files
%out = location of output files
%temp = Temperature of simulation
%pep = name of peptide
%data = input trajectory with uneven time steps
%bins used to divide the interval of the time steps
%NN = number of transitions you want gillespie traj to have
%OUTPUTS ARE:  trajectory in equal time steps

function f = tDisc(in,out,data,bins,NN,temp,pep);

P = pwd;
%Input Folder Path
myFolder = sprintf('%s\\%s', P,in);%
% Check to make sure that folder actually exists.
foldCheck(myFolder);

%Output Folder Path
OutFolder = sprintf('%s\\%s', P,out);%
foldCheck(OutFolder);

fin1 = fullfile(myFolder, data); %read in traj w/uneven tstep
data = load(fin1, '-ascii');

t = data(:,1); %time data converted to seconds
s = data(:,2); %states in second column

tim=linspace(0,t(end),bins*NN); % We discretize the time
ind=1; %index starting at 1 

state = zeros(length(tim),1);

for i=1:bins*NN
   while tim(i) > t(ind)
      ind=ind+1;
   end
   state(i)=s(ind);
end

const_t = [tim' state]; %save out time and state
fout1 = fullfile(OutFolder,sprintf('%d%sGILL_TRAJ.txt',temp,pep));
dlmwrite(fout1, const_t , 'delimiter','\t', 'precision','%.5f');

f = const_t;
end


