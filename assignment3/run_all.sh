#!/bin/tcsh
g++ -std=c++11 serial/main.cpp serial/serial_smith_waterman.cpp -o serial/serial_smith_waterman
nvcc -std=c++11 -arch=compute_52 -code=sm_52 main.cu cuda_smith_waterman_skeleton.cu -o cuda_smith_waterman
foreach i (sample.in 4k.in input1.txt input2.txt input3.txt input4.txt input5.txt input6.txt) 
  echo "> serial on datasets/$i" 
  ./serial/serial_smith_waterman "datasets/$i"
  foreach b (4 8 16)
    foreach t (32 512 1024)
      echo "> cuda on datasets/$i with $b blocks/grid and $t threads/block"
      ./cuda_smith_waterman "datasets/$i" "$b" "$t"
    end
  end
end
