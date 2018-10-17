close all


%create data set of distances (r) of each pt (x,y) to 
%the centre of each state cluster

f_in='510_dcsa_pca_12m.txt';
PCA_data=load(f_in,'-ascii') ; 

%[values, centers] = hist3([PCA_data(:,2),PCA_data(:,1)],[100,100]);
%imagesc(centers{:},values)
%colorbar
%axis square
%axis xy
plot(PCA_data(:,1), PCA_data(:,2), 'b.')


%510 ddcsa

C10 = [ 0.3, 0.914];     %possibly corelates with C10 at 450 ddcsa
C11 = [ 1.04, -0.206];   %sim position to C11 at 450 ddcsa
C12 = [-0.725, -0.302];  % possibly correlates to s4 in 450 dcsa

R = 0.15;
ns = 3;

hold on
t = linspace(0,2*pi);
plot(C10(:,1)+ R*cos(t), C10(:,2)+ R*sin(t),'g-')

axis square
axis xy

r1 = hypot(PCA_data(:,1)-C10(:,1), PCA_data(:,2)-C10(:,2));

plot(PCA_data(r1<=R,1),PCA_data(r1<=R,2),'go')

s1 = [PCA_data(r1<=R,1), PCA_data(r1 <=R,2)];

hold on 

plot(C11(:,1)+ R*cos(t), C11(:,2)+ R*sin(t),'r-')
r2 = hypot(PCA_data(:,1)-C11(:,1), PCA_data(:,2)-C11(:,2));

plot(PCA_data(r2<=R,1),PCA_data(r2<=R,2),'ro')

s2 = [PCA_data(r2<=R,1), PCA_data(r2<=R,2)];

axis square
axis xy

hold on 

plot(C12(:,1)+ R*cos(t), C12(:,2)+ R*sin(t),'k-')
r3 = hypot(PCA_data(:,1)-C12(:,1), PCA_data(:,2)-C12(:,2));

plot(PCA_data(r3<=R,1),PCA_data(r3<=R,2),'ko')

axis square
axis xy

s3 = [PCA_data(r3<=R,1), PCA_data(r3<=R,2)];

r_values = [r1, r2 ,r3];
%plot(datatest(:,1),datatest(:,2),'^')
dlmwrite('dcsa_510_rvals.txt', r_values, 'delimiter','\t') ; %write out to txt file


%PART 2

%Count Matrix Code for PCA data

%load in file with r vals from each pca pt(x,y) to the cluster centre


%
%create3 x 3 count matrix test
%this will count states only in the defined region before TBA applied
CM_test=zeros(ns,ns);

%test matrix which can id the time stamp at which transition states appear
%not sure if need this, have kept for now.
S_time=(1:length(r_values)); 

%state matrix which will use in step 3 to apply TBA, the values will be from 0-3
%this will be used later to graph the no. of states vs frames to find
%lifetimes of each state. 

State_TBA=(1:length(r_values)); %state array

%LB method uses lag time, have defined the lagtime here, set to 1 for now.
%define lagtime in n steps which is n*500ps (time between each step)
lag=1;

iter=0; %will use to determine first state i.e. s_old for the first time in
%the for loop

%figure(1);clf;hold on %plot coords for ref

%count for each state
cs1=0;
cs2=0;
cs3=0;
%cs4=0;
cstemp=0;





for i=1:lag:length(r_values)
    %for j=1:lag:length(PCA_data)

    %x=PCA_data(i,1);
    %y=PCA_data(i,2);
    
    
    %plot(x,y,'o')
    
    
    % Test for state 1 = s1
    %cartesian coords used ideally want to update to polar and use Radius
    if r_values(i,1) <= R

        snew = 1;
        cs1 = cs1 + 1;
        iter = iter + 1;

        State_TBA(i)=snew; 

        
        if iter==1
            sold = snew;
        else  
            CM_test(sold,snew) = CM_test(sold,snew) + 1;
            sold = snew;
        end

    % Test for state 2 = s2   
    elseif r_values(i,2) <= R

        snew = 2;
        cs2 = cs2 + 1;
        iter = iter + 1;
        State_TBA(i)=snew; 
        
        if iter==1
            sold = snew;
        else  
            CM_test(sold,snew) = CM_test(sold,snew) + 1;
            sold = snew;
        end  
      
    % Test for state 3 = s22   
    elseif r_values(i,3) <= R

        snew = 3;
        cs3 = cs3 + 1; 
        iter = iter + 1;
        State_TBA(i)=snew;
        
        if iter==1
            sold = snew;
        else  
            CM_test(sold,snew) = CM_test(sold,snew) + 1;
            sold = snew;
        end
        
        
        
    else
        snew= 0; %temp state
        S_time(i) = 0;
        cstemp= cstemp+1;
        State_TBA(i)=snew;
        
    end
   
        
    
end


dlmwrite('dcsa_510_state_tba.txt', State_TBA, 'delimiter','\t') ; %write out to txt file

%% PART 3
%TBA, Analysis & RATE MATRIX CREATION
state_data= State_TBA ; 
%pull out frames relating to states in non tba matrix

frame=(1:length(state_data));

State_Frames=[frame', state_data'];

%find indices to elements in first column of VLM_state that satisfy:

ind_s1= State_Frames(:,2)==1;
ind_s2= State_Frames(:,2)==2;
ind_s3= State_Frames(:,2)==3;
%ind_s4= State_Frames(:,2)==4;
%use the logical indices to index into VLM_state to return required sub
%matrices

f1= State_Frames(ind_s1,:);
f2= State_Frames(ind_s2,:);
f3= State_Frames(ind_s3,:);
%f4= State_Frames(ind_s4,:);

%write out states to txt file

dlmwrite('dcsa_510K_S1_Frames_12.txt', f1, 'delimiter','\t') ; 
dlmwrite('dcsa_510K_S2_Frames_12.txt', f2, 'delimiter','\t') ; 
dlmwrite('dcsa_510K_S3_Frames_12.txt', f3, 'delimiter','\t') ; 
%dlmwrite('dcsa_510K_S4_Frames_13.txt', f4, 'delimiter','\t') ; %write out states to txt file

%%

%creating matrix with 3 columns of data for my ref
frame=1:length(state_data);
TBA_M=[frame', S_time', state_data'];

%run length encoding 

ii= [ find(state_data(1:end-1) ~= state_data(2:end)) length(state_data) ];
len=diff([0 ii]);
val=state_data(ii);

VLM=[val', len'];



%% test may take out
%may take out - updates the time frame which is temp state to time frame of
%defined state

lastvalue=S_time(1);
z=1;
    
for x = 1:length(S_time)
    
   if S_time(1)==0
      z=z+1;
      S_time(1)=S_time(z);
   end
       
   if (S_time(x)==0) %&& (lastvalue==nextvalue)
       S_time(x)=lastvalue;
   else
       lastvalue=S_time(x);
       %nextvalue=S_time(x+2);
   end
    
end


%%

%TBA Analysis - find the temporary states and assign them to the state they
%originated from

%want to ammend this to divide the temp states in half and assign first half
%to origin and second half to next state - so far only can update to prev
%state, and the last state before the transition
%need to apply to intermediate trans states

%If the first state is transition i.e. snew=0, count forward untill find nonzero state
%update the value to the nonzero state 

State=state_data;

v=1;

if State(1)==0 

    while State(v) == 0 
         v=v+1;
         
    end
end

State(1)=State(v);


%%not necessseary check

z_idx=find(not(State));

for elem =z_idx
    
    zi= [ find(elem(1:end-1) ~= elem(2:end)) length(elem) ];
    zlen=diff([0 zi]);
    %disp(zlen)
end

%%

nz_idx=find(State); %returns index of non zero vals, not needed

%TBA

last_state=State(1);
count=0;


for y = 2:(length(State)-1)
    %for z = 2:length(State)
    
        if (State(y)==0) && State(y-1)==0
            
            State(y)=last_state;  
      
            
        elseif (State(y)==0) && (State(y+1) ~=0)
            newstate=State(y+1);

            State(y)=newstate;
            
            count=count+1;
            
        elseif (State(y)==0) && (State(y+1)==0)
            
            count=count+1;
       
            State(y)=last_state;
            
       % elseif (State(y)==0) && (State(y+1)==0) && (State(y-1)==0)
                 
        else
            last_state=State(y);
            
        end
   
   %end
    
end

%compare TBA and no TBA state arrays for reference 
%can easily see where 0 have been updated 
comp=[State_TBA' , State'];

dlmwrite('dcsa_510_TBA_check.txt', comp, 'delimiter','\t') ; %write out states to txt file
dlmwrite('dcsa_510_TBA_states.txt', State', 'delimiter', '\t');

%%run length encoding/decoding to find the lifetime of each state
%find avg number of frames per state 

num= [ find(State(1:end-1) ~= State(2:end)) length(State) ];
len1=diff([0 num]);
val1=State(num);

VLM_state=[val1', len1']; 
dlmwrite('dcsa_510_Len_Vals_check.txt', VLM_state, 'delimiter', '\t'); %write out to txt file

%find indices to elements in first column of VLM_state that satisfy:

ind1= VLM_state(:,1)==1;
ind2= VLM_state(:,1)==2;
ind3= VLM_state(:,1)==3;
%ind4= VLM_state(:,1)==4;

%use the logical indices to index into VLM_state to return required sub
%matrices

S1= VLM_state(ind1,:);
S2= VLM_state(ind2,:);
S3= VLM_state(ind3,:);
%S4= VLM_state(ind4,:);

%sum frames for each state 

S1_mean=mean(S1);
S2_mean=mean(S2);
S3_mean=mean(S3);
%S4_mean=mean(S4);

%call out second column in each mean this is mean num frames each state
%lived for i.e. mean lifetime

T1=S1_mean(2);
T2=S2_mean(2);
T3=S3_mean(2);
%T4=S4_mean(2);

Life_T=[T1 T2 T3];

%%
%Now that TBA has been applied create Count Matrix

%create 3 x 3 count matrix
CM=zeros(ns,ns);


%define lagtime in n steps which is n*500ps (time between each step)
lag_t=1;

iter_t = 0;

%figure(1);clf;hold on

%count for each state
cs1_n=0;
cs2_n=0;
cs3_n=0;
%cs4_n=0;

for w=1:lag_t:length(State)

 
    % Test for state 1 = s11
    if State(w) == 1
        s_new = 1;
        cs1_n = cs1_n + 1;
        iter_t = iter_t + 1;
        
        if iter_t==1
            s_old = s_new;
        else  
            CM(s_old,s_new) = CM(s_old,s_new) + 1;
            s_old = s_new;
        end

    % Test for state 2 = s12   
    elseif State(w) == 2
        s_new = 2;
        cs2_n = cs2_n + 1;
        iter_t = iter_t + 1;
        
        if iter_t == 1
            s_old = s_new;
        else  
            CM(s_old,s_new) = CM(s_old,s_new) + 1;
            s_old = s_new;
        end
        
    % Test for state 3 = s22   
    elseif State(w) == 3
        s_new = 3;
        cs3_n = cs3_n + 1; 
        iter_t = iter_t + 1;

        
        if iter_t ==1
            s_old = s_new;
        else
            CM(s_old,s_new) = CM(s_old,s_new) + 1;
            s_old = s_new;
        end
        
        
    end

   
end

%%

N=sum(CM); %total number of transitions from l to i where l goes from l to
%3

%gives transition matrix where PM(i,j)= CM(i,j)/sumC(i,l)
%here sum(C(i,l))==N(i)


PM=zeros(ns,ns);

PM(1,1)=CM(1,1)./(N(1));
PM(2,2)=CM(2,2)./(N(2));
PM(3,3)=CM(3,3)./(N(3));
%PM(4,4)=CM(4,4)./(N(4));

PM(1,2)=CM(1,2)./(N(1));
PM(1,3)=CM(1,3)./(N(1));
%PM(1,4)=CM(1,4)./(N(1));

PM(2,1)=CM(2,1)./(N(2));
PM(2,3)=CM(2,3)./(N(2));
%PM(2,4)=CM(2,4)./(N(2));

PM(3,1)=CM(3,1)./(N(3));
PM(3,2)=CM(3,2)./(N(3));
%PM(3,4)=CM(3,4)./(N(3));

%PM(4,1)=CM(4,1)./(N(4));
%PM(4,2)=CM(4,2)./(N(4));
%PM(4,3)=CM(4,3)./(N(4));

dlmwrite('dcsa_510_prob_matrix.txt', PM, 'delimiter', '\t'); %write out to txt file


%create rate matrix K
K=zeros(ns,ns);

%K(i,i) = -1/(T(i)) = -sum(K(i,j)
K(1,1) = -1/(Life_T(1));
K(2,2) = -1/(Life_T(2));
K(3,3) = -1/(Life_T(3));
%K(4,4) = -1/(Life_T(4));

%K(i,j) = BP(i,j)/T(i)
%branching probability of going from state i to j
%is N(i,j)sym/(sum of N(l,i)) where N is number of transitions
%the below gives output of a mathematically correct rate matrix or
%infinetismal generator, however I have not divided by the lifetime..need
%to investigate more

K(1,2)=CM(1,2)./(N(1));
K(1,3)=CM(1,3)./(N(1));
K(2,1)=CM(2,1)./(N(2));
K(2,3)=CM(2,3)./(N(2));
K(3,1)=CM(3,1)./(N(3));
K(3,2)=CM(3,2)./(N(3));


dlmwrite('510_dcsa_rate_matrix.txt', K, 'delimiter', '\t'); %write out to txt file
dlmwrite('510_dcsa_count_matrix.txt', CM_test, 'delimiter', '\t'); %write out to txt file
dlmwrite('510_dcsa_count_matrix_TBA.txt', CM, 'delimiter', '\t'); %write out to txt file

%%
%find stationary distribution of K which is the equilibrium pop
%singular value decompositing to get peq that satisifies K*peq=0
%this is not working, need to get non trivial answer here, matlab returning
% all zeros. Can work out by hand but not ideal

%[U S V] = svd(inv(K));
%peq=V(:,end);

%%PART 4

states=State' ;  

%create time vector 

time = linspace(0,length(states)*.0005,length(states));
t = time';

st = [t , states];
dlmwrite('dcsa_510K_TBA_traj.txt', st, 'delimiter','\t') ;

figure() 

plot(t(:,1), states(:,1), 'b-')

xlabel('time (\mu s) ')
ylabel('state')

xlim([0 3])
ylim([0 3.5])
xticks([0 1 2 3])
yticks([ 1 2 3 ])
