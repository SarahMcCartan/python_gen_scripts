SLIDING WINDOW COUNT MATRIX

INPUT = DISCRETE TBA TRAJECTORY FOR BEST R VALUE 

slide frames looking at by specific lag, lag = 10 initially
slide per Y number of frames, set Y = 1 initially
data saved out every 0.5ns or 1ns in 2 cases
del = 0.5 or 1 ns

S =state data:

n=0: S(0) -> S(lag)
n=1: S(del)_> S(del+lag)
n=2: S(2*del) -> S(2*del + lag)
etc

OUTPUT = Count Matrix, only interested in the transitions
Number of counts (transitions and non transitions) = length(State) - lag


CM(1,2) = Number of Transitions from S2 to S1
CM(1,1) = Number of Counts that S1 stayed in S1


Next Step Symmetrise COUNT MATRIX, THEN make estimate of the TM to get estimate of  rel time 
from the first non zero eval

