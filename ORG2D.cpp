#include<vector>
#include<cmath>
#include<cstdlib>
#include<cstdio>
#include<set>
#include<ctime>
#include<cassert>
#include<algorithm>
#include<numeric>
#include <chrono> 
#include <iostream>   
#include <fstream>   
#include<functional>
#include "gurobi_c++.h"
#include"vec2.h"
#include<omp.h>

#define GS 64 // Grid Scale: GS x GS

vec2l farthest_clustering(const vec2l &pts, int K, double& radius){
    // Farthest clustering for pts in O(number of samples * number of centers)
    // Return value: center lists
    // radius: resulting radius

	int n = pts.size();
	if(K>=n) {radius= 0; return pts;}
	vector<double> dist_to_nearest_point(n, std::numeric_limits<double>::max());
    vector<bool> chosen(n, false);
    vec2l res;
    int max_dist_id = 0;
    double max_val;
    while(K--){
        chosen[max_dist_id] = 0;
        auto cen = pts[max_dist_id];
        res.push_back(cen);
        max_val = 0.0;
        for(int i=0;i<n;i++){
            if(!chosen[i]){
                dist_to_nearest_point [i] = min(dist_to_nearest_point[i], dis(cen, pts[i]));
                if(dist_to_nearest_point[i] > max_val){
                    max_val = dist_to_nearest_point[i];
                    max_dist_id = i;
                }
            }
        }
    }
    radius = max_val;
	return res;
}

vec2l iterative_centering(const vec2l &pts, int K, double& radius){
    // use iterative centering
    // step1, sample K random points
    // step2, cluster will 1-nearest neighbor
    // step3, apply 1-center for each cluster, and go to step2 unless center does not update anymore
    // resample several times
    vec2l pts_local = pts;
    shuffle(begin(pts_local), end(pts_local));
    vec2l centers;
    for(int i=0;i<K;i++) centers.push_back(pts_local[i]);
    while(1){
        vector<vector<int>> grouped_pts(K);
        // voronoi
        
    }
}

vec2l ilp(const vec2l &pts, const polygon &p, int K, double &radius){
    GRBEnv env = GRBEnv(true);
    env.set(GRB_IntParam_OutputFlag, false);
    env.start();
    int m = pts.size();
    vector<int> xco, yco;
    vec2l res;

    for(auto &v: pts){
        xco.push_back((int)(v.x * GS));
        yco.push_back((int)(v.y * GS));
    }
    double ma = sqrt(2) / 2, mi = p.circumference / 2 / m;
    double thres = 1.0 / GS / 4;
    
    GRBVar (*ys)[GS] = new GRBVar[GS][GS];
    do{
        GRBModel model = GRBModel (env);
        double mid = (ma + mi) / 2;
        for(int i=0; i < GS; i++){
            for(int j=0; j < GS; j++){
                ys[i][j] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY, "y_" + to_string(i) + "_" + to_string(j));
            }
        }
        for(int i=0; i < m; i++){
            int x = xco[i];
            int y = yco[i];
	        int xmin = max(0, int(floor(xco[i] - mid * GS)));
            int ymin = max(0, int(floor(yco[i] - mid * GS)));
            int xmax = min(GS-1, int(ceil(xco[i] + mid * GS)));
            int ymax = min(GS-1, int(ceil(yco[i] + mid * GS)));
            GRBLinExpr expr = 0;
            for (int j = xmin; j <= xmax; j++) {
                for (int k = ymin; k <= ymax; k++) {
                    if ((j - x)*(j - x) + (k - y) * (k - y) <= GS * GS * mid * mid) {
                        expr += ys[j][k];
                    }
                }
            }
            model.addConstr(expr >= 1, "cover_" + to_string(i));
        }
        GRBLinExpr expr = -K;
        for (int i = 0; i < GS; i++) {
            for (int j = 0; j < GS; j++) {
                expr += ys[i][j];
            }
        }
        model.addConstr(expr <= 0, "centers_constr");
        model.optimize();
        if (model.get(GRB_IntAttr_SolCount) > 0) {
            ma = mid;
            radius = ma;
            res.clear();
            for(int i=0; i<GS; i++){
                for(int j=0; j<GS; j++){
                    if(ys[i][j].get(GRB_DoubleAttr_X) == 1.0){
                        res.emplace_back((i+.5)/GS, (j+.5)/GS);
                    }
                }
            }
        }
        else {
            mi = mid;
        }
    }while(ma - mi > thres);
        
    delete[] ys;
    return res;
}

vec2l sample2d(polygon& p, int gs){
    // put the polygon into a GS x GS grid in 1 x 1 square
    vec2l res;
    double cl = 1.0/gs;
	int n = p.pts.size();
    for(double ch = cl/2; ch < 1; ch += cl){
        vector<double> seq;
        for(int j=0;j<n;j++){
			if(p.pts[j].y > ch && p.pts[(j+1)%n].y<=ch || p.pts[j].y<ch && p.pts[(j+1)%n].y>=ch ){
				seq.push_back(p.pts[j].x + (p.pts[(j+1)%n].x - p.pts[j].x) * (ch - p.pts[j].y) /  (p.pts[(j+1)%n].y - p.pts[j].y) );
			}
        }
        sort(begin(seq), end(seq));
        while(seq.size() > 0){
            double r = seq.back();
            seq.pop_back();
            double l = seq.back();
            while(r >= l){
                res.emplace_back(floor(r/cl) * cl + cl / 2, ch);
                r -= cl;
            }
            seq.pop_back();
        }
    }
    return res;
}

int main(int argc, char*argv[]){
    if(argc < 2){
        printf("Format: ./org2d algo_number < poroblem.txt\n\t"
        "1: Farthest clustering\n\t"
        "2: Integer Progamming Solver\n");
        return 0;
    }
    // Number of vertex in the polygon, number of samples (does't count for ORG2D), number of centers(robots)
    int n, m, k; 

    polygon p;
    double x, y;
    double x_max = numeric_limits<double>::lowest(), x_min = numeric_limits<double>::max();
    double y_max = numeric_limits<double>::lowest(), y_min = numeric_limits<double>::max();

    scanf("%d%d%d", &n, &m, &k); 
    for(int i=0;i<n;i++){
        scanf("%lf%lf", &x, &y);
        x_min = min(x_min, x); x_max = max(x_max, x);
        y_min = min(y_min, y); y_max = max(y_max, y);
        p.pts.emplace_back(x, y);
    }
    
    // polygon is scaled to be contained in a 1x1 square
    double scale = max(x_max - x_min, y_max - y_min);
    for(auto &v: p.pts){
        v.x = (v.x - x_min) / scale;
        v.y = (v.y - y_min) / scale;
    }

    vec2l centers;
    double radius;
    
    if(argv[1][0] == '1'){
        auto sam_pts = sample2d(p, 1000);
        centers = farthest_clustering(sam_pts, k, radius);
    }else if(argv[1][0] == '2'){
        auto sam_pts = sample2d(p, GS);
        centers = ilp(sam_pts, p, k, radius);
    }
    // rescale back
    radius *= scale;
    for(auto &v: centers){
        v.x = v.x * scale + x_min;
        v.y = v.y * scale + y_min;
    }
    printf("%.4lf\n", radius);
    for(auto v: centers){
        printf("%.4lf %.4lf\n", v.x, v.y);
    }
    return 0;
}
