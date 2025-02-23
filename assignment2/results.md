|          | sample.in   | 1k.in     | 20k.in  | 20k2k.in | input1.txt  | input2.txt  | input3.txt |
| -------- | ----------- | --------- | ------- | -------- | ----------- | ----------- | ------- |
| *Answer* | 13          | 662       | 14333   | 1565     | 22          | 9           | 35246   |
| serial   | 8.641e-06 | 0.0268735 | 10.667 | 1.08233 | 1.7736e-05 | 9.788e-06 | 6.97727 |
| n=1      | 7.511e-05 | 0.0316708 | 12.2597 | 1.26028 | 9.1409e-05 | 7.6577e-05 | 8.28596 |
| n=2      | 9.1173e-05 | 0.0169174 | 6.29639 | 0.667237 | 0.000114075 | 9.6628e-05 | 4.28979 |
| n=4      | 0.000424237 | 0.0154124 | 5.2209 | 0.590083 | 0.00048403 | 0.000220468 | 3.57744 |
| n=7      | 0.000412461 | 0.0195966 | 6.23499 | 0.719763 | 0.000809093 | 0.00074339 | 4.26571 |
| n=8      | 0.00053384 | 0.0192869 | 5.76089 | 0.684382 | 0.00094113 | 0.000835815 | 3.96675 |



Lesson

- Use local variables as much as possible. Use global arrays as restrictive as possible
- Memory access pattern is the speedup bottleneck
- Critical sections can lurk around

From solution

- Barrier-only is the way
- It is completely okay to process element-by-element, as long as no **file-global variables** are involved, it will probably be fast enough
- As a middle ground, declare main-local variables (pointers) and pass them into thread_word (as an optimisation this can be done as a struct to avoid further per-thread derviations and passing many arguments into thread_work)