### 2D coverage

The problem asks to find **k** static locations on the plane to minimize the coverage radius to guard a polygon    Optimal Region Guarding (ORG2D), Optimal Perimeter Guarding(OPG2D) with 2D circles (robots with 2D range sensors).

The problem is hard, even to approximate within a factor of 1.152 even for a simple polygon.

This repo contains algorithms described in paper https://arxiv.org/pdf/2002.08477.pdf in R:SS 2020



To run it,

1. install GuroBi, and update the GuroBi installation location in makefile
2. make all
3. run `./opg2d algo_number < problem.txt > ans.txt` for OPG2D or `./org2d algo_number < problem.txt > ans.txt` for ORG2D
   - `problem.txt` should contain the number of vertex, the number of samples, the number of centers (robots), and then the polygon coordinates in clockwise or counter-clockwise order.
4. visualize the solution using `python draw_algo problem.txt ans.txt`