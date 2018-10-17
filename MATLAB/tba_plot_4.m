f_in = 'dcsa_510k_TBA_traj.txt';

states=load(f_in,'-ascii') ;  

%create time vector 

time = linspace(0,length(states)*.001,length(states));
t = time';

st = [t , states];
%dlmwrite('CsA_310K_TBA_traj.txt', st, 'delimiter','\t') ; 

plot(t(:,1), states(:,2), 'b-')

xlabel('time (\mu s) ')
ylabel('state')

xlim([0 3])
ylim([8.5 11.5])
xticks([1 2 3])
yticks([ 9 10 11 ])
