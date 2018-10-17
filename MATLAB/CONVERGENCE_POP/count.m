%% FUNCTION TO PERFORM SLIDING WINDOW COUNTING FOR AN ARRAY OF LAG TIMES
%in = location of input files
%out = location of output files
%temp = Temperature of simulation
%pep = name of peptide
%ns = number of states
%data = input TBA trajectory
%a:b:c = start:interval:end to give array of lag times
%frame is number of frames to apply sliding window to 
%sd = col where state data is
%check = 1 save out check index, = 0 doesnt
%OUTPUTS ARE:   Count_Matrix_TBA = Matrix of No. of Transitions per lag value
%               Check_Index = check to ensure sliding window applied correctly  
%               lagVals = all the lag values saved out

%%
function CM = count(in, out,temp, pep,ns, data,a,b,c, frame,sd,check)
P = pwd; %find current directory

%INPUT FOLDER PATH
myFolder = sprintf('%s\\%s', P,in);
% Check to make sure that folder actually exists.  Warn user if it doesn't.
foldCheck(myFolder);

%OUTPUT FOLDER PATH
OutFolder = sprintf('%s\\%s', P,out);
foldCheck(OutFolder);
    
    fin = load(sprintf('%s\\%s', myFolder,data),'-ascii');

    State = fin(:,sd); %read in col of state data
    iter = 0;
    lenMat = floor(c/b); %find number of lag values 
    lagVals = zeros(lenMat , 1); %empty col vector
    C=zeros(ns,ns); %empty matrix

    for lag = a:b:c

        iter = iter+1;
        lagVals(iter) =lag; %save out lag values

        %sliding window
        % TBA(1) -> TBA(lag), TBA(2) -> TBA(2+lag)
        % generally: TBA((n+1)frame)) -> TBA((n+1)*frame + lag))

        w1 = zeros((length(State)-lag-frame+1),1);
        w2 = zeros((length(State)-lag-frame+1),1);

        for n = 0:(length(State)-lag-frame)

            for j=1:2

                if(j==1)
                    w=n+1;
                    iter_t = 1;
                    w1(n+1)=w;
                else
                    w=n+1+lag;
                    iter_t = 2;
                    w2(n+1)=w;
                end

                for st = 1:ns

                    if State(w) == st
                        s_new = st;

                        if iter_t==1 %condition for first state visited
                            s_old = s_new;
                        else
                            C(s_old,s_new) = C(s_old,s_new) + 1;
                            s_old = s_new; %update state to count tranisition
                        end

                    end

                end


            end
        end
        
        CM = C;
        fout1 = fullfile(OutFolder, sprintf('%d%s_CM%06d.txt',temp,pep,lag));
        dlmwrite(fout1, CM, 'delimiter', '\t', 'precision','%.5f');
        
        if check ==1
            
            fout2 = fullfile(OutFolder, sprintf('%d%sCheckIndex%06d.txt',temp,pep,lag));
            dlmwrite(fout2, [w1 w2], 'delimiter', '\t', 'precision','%.5f');
        end

        C=zeros(ns,ns); %reset the CM for next loop

    end

    fout3 = fullfile(OutFolder, sprintf('%d%s_lagVals.txt',temp,pep));
    dlmwrite(fout3, lagVals, 'delimiter', '\t');
    
   
    end
