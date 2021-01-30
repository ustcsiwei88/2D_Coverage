### 2D coverage

The problem asks to find **k** static locations on the plane to minimize the coverage radius,
to guard 
1. a polygon boundary; Optimal Perimeter Guarding (OPG2D) with 2D circles (robots with 2D range sensors).
2. a polygonal region; Optimal Region Guarding with 2D circles (ORG2D), 

The problem is NP-hard, even to approximate it within a factor of 1.152 for a simple polygon.

This repo contains algorithms described in paper https://arxiv.org/pdf/2002.08477.pdf published in R:SS 2020

To run and visualize the result of the algorithms, please take the following steps

1. install Gurobi, and update Gurobi installation location in Makefile
2. run `make all`
3. run `./opg2d algo_number < problem.txt > ans.txt` for OPG2D or `./org2d algo_number < problem.txt > ans.txt` for ORG2D
   - `problem.txt` should contain the number of vertex, the number of samples, the number of centers (robots), and then the polygon coordinates in clockwise or counter-clockwise order.
4. visualize the solution using `python draw_algo problem.txt ans.txt`
