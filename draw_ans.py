import sys
import matplotlib.pyplot as plt

if __name__ == '__main__':
    if len(sys.argv) < 3:
        print("Format:\n\t python draw_ans problem.txt ans.txt")
        exit(0)
    problem = open(sys.argv[1])
    ans = open(sys.argv[2])

    fig, ax = plt.subplots()

    n = int(problem.readline()) # polygon vertices number
    m = int(problem.readline()) # sample number
    k = int(problem.readline()) # centers (robots) number

    x_coor, y_coor = [], []
    for i in range(n):
        line = problem.readline().split(' ')
        x_coor.append(float(line[0]))
        y_coor.append(float(line[1]))
    x_coor.append(x_coor[0])
    y_coor.append(y_coor[0])
    plt.plot(x_coor, y_coor, 'k-')
    
    radius = float(ans.readline())
    cen_x, cen_y = [], []
    for i in range(k):
        line = ans.readline().split(' ')
        cen_x.append(float(line[0]))
        cen_y.append(float(line[1]))
        circle_ =  plt.Circle((cen_x[-1], cen_y[-1]), radius, fill=False, color='r')
        ax.add_patch(circle_)

    ax.axis('off')
    ax.axis('equal')
    plt.show()
