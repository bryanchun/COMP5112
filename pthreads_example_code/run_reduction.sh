echo "reduce algorithm 0"
nvcc -o reduction0 -D REDUCE0 -arch=compute_52 -code=sm_52 reduction.cu -run
echo " "
echo "reduce algorithm 1"
nvcc -o reduction1 -D REDUCE1 -arch=compute_52 -code=sm_52 reduction.cu -run
echo " "
echo "reduce algorithm 2"
nvcc -o reduction2 -D REDUCE2 -arch=compute_52 -code=sm_52 reduction.cu -run
echo " "
echo "reduce algorithm 3"
nvcc -o reduction3 -D REDUCE3 -arch=compute_52 -code=sm_52 reduction.cu -run
