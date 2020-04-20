#!/bin/tcsh
g++ -std=c++11 serial/main.cpp serial/serial_smith_waterman.cpp -o serial/serial_smith_waterman
g++ -std=c++11 -lpthread main.cpp pthreads_smith_waterman_skeleton.cpp -o pthreads_smith_waterman
foreach i (sample.in 1k.in 20k.in 20k2k.in input1.txt input2.txt input3.txt) 
  echo "> serial on datasets/$i" 
  ./serial/serial_smith_waterman "datasets/$i"
  foreach p (1 2 4 7 8)
    echo "> pthread on datasets/$i with $p threads"
    ./pthreads_smith_waterman "datasets/$i" "$p"
  end
end
