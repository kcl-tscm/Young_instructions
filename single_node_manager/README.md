This is a file for submitting within the same node a set of 450 files (this can be changed by user), one calculation is assigned to one core. In young you have 40 cores. 
The number of files (450) will be partioned and each calculation once finished will be completed as the ./Cos par_$i disposes the output.


To use it in Young, just copy and paste the whole script into a file called job.sh
Change the \$ variables accordingly. 
use the command chmod +x job.sh for enabling running it.
