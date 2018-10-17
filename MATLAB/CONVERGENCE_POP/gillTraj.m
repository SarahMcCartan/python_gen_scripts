
%% GILLESPIE ALGORITHM TO GENERATE LONGER TIME SERIES DATA FROM RATE MATRIX
%in = location of input files
%out = location of output files
%temp = Temperature of simulation
%pep = name of peptide
%ns = number of states
%data = input rate matrix to generate traj from
%NN = number of transitions you want gillespie traj to have
%OUTPUTS ARE:  gillespie trajectory between states with uneven time step
% uneven as gillTraj finds the time that a transition will occur.

function f = gillTraj(in,out,data,ns,NN,temp,pep)
P = pwd;

%Input Folder Path
myFolder = sprintf('%s\\%s', P,in);%
% Check to make sure that folder actually exists.
foldCheck(myFolder);

%Output Folder Path
OutFolder = sprintf('%s\\%s', P,out);%
foldCheck(OutFolder);

fin1 = fullfile(myFolder, data); %read in rate matrix
K = load(fin1, '-ascii');

%Spectral decomposition (getting eigenvalues/vectors) and
%finding equilibrium probability vector')
%calculate equilibrium from spectral decomposition

[eigvec,eigval]=eig(K); % diagonalize K, eigvec stores the eigenvectors, eigval the eigenvalues
[dsorted,index]=sort(diag(eigval),'descend');
% sort the eigenvalues. dsorted stores the eigenvalues, index the corresponding indices

% Plot of sorted eigenvalues
figure
hold on
x=linspace(1,ns,ns);
plot(x,dsorted,'o')
ylabel('Eigenvalue','FontSize',18)
% sorted eigenvalues:
ind=index(1); %index corresponding to 0 eigenvalue

Peq=eigvec(:,ind)/sum(eigvec(:,ind));
%equilibrium probability corresponds to 0 eigenvalue.
%Based on the equilibrium probability we can also obtain the energy.

%Finding second right eigenvector

% splitting and eigvec
[eigvec,eigval]=eig(K'); % diagonalize K, eigvec stores the right eigenvectors
[dsorted,index]=sort(diag(eigval),'descend'); % sort the eigenvalues.

slowest_relrate=-dsorted(2);
slow_vec=eigvec(:,index(2));


%% Run Gillespie to obtain trajectories
for i=1:ns
    for j=1:ns
        if ( i ~= j ) % i =\= j
            p(j,i)=K(j,i)/(-K(i,i));
        end
    end
end

for j=1:ns
    pp(1,j)=0;
    for i=1:ns
        pp(i+1,j)=sum(p(1:i,j)); % The pp matrix stores cumulative transition probabilities
    end
end

s=1; %starting state
time=0;
tcum=zeros(ns,1); % creates a vector of ns elements

for k=1:NN
    tadd=(1/K(s,s))*log(rand); % tadd is the survival time the system stays in state s.
    %It is calculated based on a single exponential decay for the survival time.
    time=time+tadd; %time=total time
    t_traj(k)=time; %t_traj matrix saves information about the survival times
    s_traj(k)=s;    %s_traj matrix saves information about the current state
    ss=find(histc(rand,pp(:,s))); % finds the new state
    tcum(s)=tcum(s)+tadd; % total time that the system spends in state s
    s=ss; %the new state is ss
end

states_time = [ t_traj' s_traj'];
fout1 = fullfile(OutFolder, sprintf('%d%s_gillTraj_rough.txt',temp,pep));
dlmwrite(fout1, states_time, 'delimiter', '\t');

%average time spent in each state / total time = measured euqilibrium probability (p_eq)
fprintf('Comparing the true eigenvector to the trajectories estimate')
measured_Peq=tcum/sum(tcum);
display(Peq)
display(measured_Peq)

f = states_time;
end


