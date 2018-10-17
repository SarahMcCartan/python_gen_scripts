clear variables


P = pwd;
in = 'functionTest\';
fname = 'KMRND3.txt';
ns =3 ;

% Specify the folder where the files live.
myFolder = sprintf('%s\\%s\\', P, in);
foldCheck(myFolder); %check folder exists

fin = fullfile(myFolder, fname); %read in lag time step file
K = load(fin, '-ascii');

fprintf ('calculate equilibrium from spectral decomposition\n')
[eigvec,eigval]=eig(K); % diagonalize K, eigvec stores the eigenvectors, eigval the eigenvalues

display(eigval)

[dsorted,index]=sort(diag(eigval),'descend'); % sort the eigenvalues. dsorted stores the eigenvalues, index the corresponding indices
% sorted eigenvalues:
ind=index(1); %index corresponding to 0 eigenvalue
eq=eigvec(:,ind)/sum(eigvec(:,ind)); % equilibrium probability corresponds to 0 eigenvalue. Based on the equilibrium probability we can also obtain the energy.

display(eq)

fprintf(' Plot of sorted eigenvalues\n')
figure
hold on
x=linspace(1,ns,ns);
plot(x,dsorted,'o')
ylabel('Eigenvalue','FontSize',18)
%% 

fprintf('Calculating equilibrium energies\n ')

kB=0.0019872041; % Boltzmann constant (kcal/mol)
temp=298; % Temperature
energy=kB*temp*(-log(eq));% calculate the energy
energy=energy-min(energy); %
% Now we plot the energies...
figure
hold on
xlabel('# State','FontSize',18)                   
ylabel('\DeltaG (kcal/mol)','FontSize',18)
bar(energy,'r')
% plot(energy,'b-o','MarkerSize',10)
% hold off

%% 

fprintf(' Finding second right eigenvector\n')

% splitting and eigvec
[eigvec,eigval]=eig(K'); % diagonalize K, eigvec stores the right eigenvectors
[dsorted,index]=sort(diag(eigval),'descend'); % sort the eigenvalues. 

display(dsorted)

slowest_relrate=-dsorted(2);
slow_vec=eigvec(:,index(2));

figure
hold on
bar(slow_vec)
xlabel('# State','FontSize',18)
ylabel('Second eigenvector','FontSize',18)
hold off
%% 
fprintf('Spectral Decomposition Option 2\n')

 %sort eigenvalues and vectors in descending order s.t. eval = 1 is first

[V, D, W] = eig(K); %V = right evectors, D = eigenvalues, W = left evectors

[D,idx] = sort(diag(D), 'descend');
disp(D)
V = V(:, idx);
disp(V)
W = W(:,idx);
disp(W)

%%
fprintf('transpose K\n')

[V2, D2, W2] = eig(K'); %V = right evectors, D = eigenvalues, W = left evectors

[D2,idx2] = sort(diag(D2), 'descend');
disp(D2)
V2 = V2(:, idx2);
disp(V2)
W2 = W2(:,idx2);
disp(W2)

