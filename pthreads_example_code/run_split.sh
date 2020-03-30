echo "run split algorithm"
nvcc -I /usr/local/cuda/samples/common/inc -o split -run -arch=compute_52 -code=sm_52 split.cu 
