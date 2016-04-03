#!/bin/bash

function runloops {
  for nants in 1 10 ; do
    for nsteps in 10 100; do
      randprob=0.4
      evaprate=0.1
      ants capitals Olympia Augusta $nants $nsteps $randprob $evaprate > o
      pathlength=`tail -1 o`
      echo  $pathlength $nants $nsteps $randprob $evaprate
    done
  done
}

function seeresults {
  sort -n $1 > sorted
  echo "Results in file named \"sorted\""
  awk 'BEGIN {print "SolutionLength nAnts nSteps randProb evapRate"}\
             {printf "%10d %9d %5d %8.2f %7.3f\n", $1, $2, $3, $4, $5}' sorted | more
}

cat <<EOF
Functions defined:
  runloops > results    #runs multiple loops and directs results to file "results"
  seeresults results    #sorts and prints results
EOF


