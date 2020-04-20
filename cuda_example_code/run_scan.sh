echo "Prefix Scan1"
nvcc -o scan_example1 -run -D SCAN1 -arch=compute_52 -code=sm_52 scan.cu
echo ""
echo "Prefix Scan2"
nvcc -o scan_example2 -run -D SCAN2 -arch=compute_52 -code=sm_52 scan.cu
