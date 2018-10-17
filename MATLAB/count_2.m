close all 
clear all

%Count Matrix Code for PCA data

%load in file with r vals from each pca pt(x,y) to the cluster centre

f_in='csa_310_rvals.txt';
PCA_data=load(f_in,'-ascii') ;  

%define Max Radius from cluster centre to be defined in that state
R = 0.15; 
%%
%create3 x 3 count matrix test
%this will count states only in the defined region before TBA applied
ns=3;
CM_test=zeros(ns,ns);

%test matrix which can id the time stamp at which transition states appear
%not sure if need this, have kept for now.
S_time=(1:length(PCA_data)); 

%state matrix which will use in step 3 to apply TBA, the values will be from 0-3
%this will be used later to graph the no. of states vs frames to find
%lifetimes of each state. 

State_TBA=(1:length(PCA_data)); %state array

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





for i=1:lag:length(PCA_data)
    %for j=1:lag:length(PCA_data)

    %x=PCA_data(i,1);
    %y=PCA_data(i,2);
    
    
    %plot(x,y,'o')
    
    
    % Test for state 1 = s1
    %cartesian coords used ideally want to update to polar and use Radius
    if PCA_data(i,1) <= R

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
    elseif PCA_data(i,2) <= R

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
    elseif PCA_data(i,3) <= R

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


dlmwrite('csa_310_state_tba.txt', State_TBA, 'delimiter','\t') ; %write out to txt file
