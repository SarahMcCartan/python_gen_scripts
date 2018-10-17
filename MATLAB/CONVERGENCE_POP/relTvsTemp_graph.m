clear variables


fin = 'dCSA_reltvsTemp.txt'; %read in lag file
CSA = load(fin, '-ascii');

y = CSA(:,1);
x = CSA(:,2);

scatter(x , y,'filled')

ylabel('Temp (K)')
xlabel('Relaxation time (s)')


%ylim([0 70])

%%
%dcsa
