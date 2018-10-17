f_in = 'csa_310_TBA_states.txt';

states=load(f_in,'-ascii') ;  

time = linspace(0,3.03,length(states));
t = time';

st = [t , states];
dlmwrite('CsA_310K_TBA_traj.txt', st, 'delimiter','\t') ; 

plot(t(:,1), states(:,1), 'b-')

xlim([0 3])
ylim([0 3.5])
xticks([0 1 2 3])
yticks([ 1 2 3 ])
