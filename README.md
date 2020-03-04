# COMP5112
HKUST COMP 5112 Spring 2020 Parallel Programming, Code Repository

### MPI

```bash
mpicc [-std=c99] -g -Wall -o $MPI_FILE $MPI_FILE.c
mpiexec [-oversubscribe] [-hostfile $HOST_FILE] -n $NUM_PROCESSORS ./$MPI_FILE
```
