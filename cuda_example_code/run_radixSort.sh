echo "run radixSort algorithm"
nvcc -I /usr/local/cuda/samples/common/inc -o radixSort -run -arch=compute_52 -code=sm_52 radixSort.cu
