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

#define GS 200 // Grid Scale: GS x GS

vec2l farthest_clustering(const vec2l &pts, int K, double& radius){
    // Farthest clustering for pts in O(nk)
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

vec2 enc_circle(vec2 p1, vec2 p2, vec2 p3, double &r){
    // Enclosing circle of the three points
    // Return: center
    // r: radius
    vec2 cen;
	double a = dis(p2, p3), b = dis(p3, p1), c = dis(p1, p2);
	r = a * b * c / sqrt((a+b+c)*(a+b-c)*(a+c-b)*(b+c-a));
	a = norm(p1), b= norm(p2), c= norm(p3);
	double d = (p1.x * (p2.y - p3.y) + p2.x * (p3.y - p1.y) + p3.x * (p1.y - p2.y)) * 2;
	cen.x = (p1.y*(c - b) + p2.y*(a - c) + p3.y*(b - a)) / d;
	cen.y = (p1.x*(b - c) + p2.x*(c - a) + p3.x*(a - b)) / d;
    return cen;
}

vec2 mini_disc(vec2l &pts, double &r){
    // minimum disc containing pts
	int n = pts.size();
    random_shuffle(begin(pts), end(pts));
	r = dis(pts[1], pts[0])/2;
	vec2 cen = (pts[0] + pts[1])/2;
	for(int i=2;i<n;i++){
		if(dis(cen, pts[i]) <= r) continue;
		// mini_disc_with_1_point: pts[i]
        r = dis(pts[0], pts[i])/2;
        cen = (pts[0] + pts[i])/2;
        for(int j=1;j<i;j++){
            if(dis(cen, pts[j]) <= r) continue;
            // mini_disc_with_2_points: pts[i] pts[j]
            r = dis(pts[i], pts[j])/2;
            cen = (pts[i] + pts[j])/2;
            for(int k=0; k<j; k++){
                if(dis(cen, pts[k]) <= r) continue;
                cen = enc_circle(pts[i], pts[j], pts[k], r);
            }
        }
	}
    return cen;
}

void compute_rmost(const vec2l &pts, vector<int>& rmost, double r){
    // rmost[i]: most number of vertices radius r circle can tile from pts[i]
    int end = 1, n = pts.size();
    vec2l seq{pts[0], pts[1]};
    double cur_r;
    for(int i=0; i<n; i++){
        mini_disc(seq, cur_r);
        while(cur_r <= r){
            end = (end+1) % n;
            if(seq.size() == n){
                // circle contains all vertices
                fill(std::begin(rmost), std::end(rmost), n);
                return ;
            }
            seq.push_back(pts[end]);
            mini_disc(seq, cur_r);
        }
        rmost[i] = (end - i + n ) % n; 
        for(auto it = begin(seq);; ++it){
            if(*it == pts[i]){seq.erase(it); break;}
        }
    }
}

vec2l contin_(const vec2l &pts, const polygon &p, int K, double &radius){
    // 1 + epsilon solution for covering perimeter with continuous coverage assumption
    // pts: sampled points on the perimeter
    // K: number of centers
    // O(n^2 log(n)) but O(n^2 log(n) / K^2) on average
    int n = pts.size(), cnt, left, pos;
    double ma = p.circumference / 2 / K;
    double mi = p.circumference / 2 / n;
    double thres = p.circumference / n / 10, r;
    vector<int> rmost(n);
    vec2 cen;
    vec2l res; // center list
    if(K==1){
        auto copy_pts = pts;
        return {mini_disc(copy_pts, radius)};
    }
    if(K >= n){
        radius = p.circumference / n / 2;
        return pts;
    }
    do{
        start_:
        double mid = (ma + mi)/2;
        compute_rmost(pts, rmost, mid);
        // Try different tiling locations
        for(int i=0; i<n; i++){
            int cnt = 0, left = K;
            int pos = i;
            while(left -- ){
                if((cnt += rmost[pos]) >= n){
                    ma = mid;
                    goto start_;
                }
                (pos += rmost[pos]) %= n;
            }
        }
        mi = mid;
    }while(ma - mi > thres);
    // Generating center list
    compute_rmost(pts, rmost, ma);
    for(int i=0; i<n; i++){
        int cnt = 0, left = K;
        int pos = i;
        while(left -- ){
            if((cnt += rmost[pos]) >= n){
                radius = 0;
                pos = i;
                cnt = 0;
                do{
                    vec2l seq;
                    int n_pts = rmost[pos];
                    for(int j = 0; j < n_pts; j++){
                        seq.push_back(pts[pos]);
                        if(++pos == n) pos = 0;
                        if(pos == i) break;
                    }
                    res.push_back(mini_disc(seq, r));
                    radius = max(radius, r);
                }while(pos != i);
                radius += p.circumference / 2 / n;
                return res;
            }
            (pos += rmost[pos]) %= n;
        }
    }
    return {};
}

vec2l ilp(const vec2l &pts, const polygon &p, int K, double &radius){
    GRBEnv env = GRBEnv(true);
    env.set(GRB_IntParam_OutputFlag, false);
    env.start();
    int m = pts.size();
    vector<int> xco, yco;
    vec2l res;
    for(auto &v: pts){
        xco.push_back(int(v.x * GS));
        yco.push_back(int(v.y * GS));
    }
    double ma = p.circumference / K / 2 + 1.0 / GS, mi = p.circumference / 2 / m;
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
            GRBLinExpr expr ;
            for (int j = xmin; j <= xmax; j++) {
                for (int k = ymin; k <= ymax; k++) {
                    if ((j - x) * (j - x) + (k - y) * (k - y) <= GS * GS * mid * mid) {
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

vec2l sample(polygon& p, int N){
    // sample N points from the boundary of the polygon
    vec2l res;
    double d = .0, td = .0;
    int n = p.pts.size(), pos = 0;
    if(p.dis.size() == 0){
        for(int i=1;i<n;i++){
			p.dis.push_back(dis(p.pts[i], p.pts[i-1]));
		}
        p.dis.push_back(dis(p.pts[0], p.pts.back()));
        p.circumference = accumulate(begin(p.dis), end(p.dis), 0.0);
    }
    d = p.circumference / N;
    res.push_back(p.pts[0]);
    while(--N){
        td += d;
        while(td > p.dis[pos]){
            td -= p.dis[pos];
            (pos+=1) %= n;
        }
        res.emplace_back(p.pts[pos] + (p.pts[(pos+1)%n] - p.pts[pos]) * (td / p.dis[pos]));
    }
    return res;
}

int main(int argc, char*argv[]){
    if(argc < 2){
        printf("Format: ./opg2d algo_number < poroblem.txt\n\t"
        "1: Farthest clustering\n\t"
        "2: Coverage with continuous coverage limitation\n\t"
        "3: Integer Progamming Solver\n");
        return 0;
    }
    int n, m, k; // Number of vertex in the polygon, number of samples, number of centers(robots)
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
    vec2l centers;
    double radius;
    if(argv[1][0] == '1'){
        auto sam_pts = sample(p, n);
        centers = farthest_clustering(sam_pts, k, radius);
    }else if(argv[1][0] == '2'){
        auto sam_pts = sample(p, n);
        centers = contin_(sam_pts, p, k, radius);
    }else{
        // polygon is scaled to be contained in a 1x1 square
        double scale = max(x_max - x_min, y_max - y_min);
        for(auto &v: p.pts){
            v.x = (v.x - x_min) / scale;
            v.y = (v.y - y_min) / scale;
        }
        // each length O(1/GS) segment along the perimeter needs a sample, 
        // so set to approximately 2 * Grid Scale
        auto sam_pts = sample(p, GS * 2); // sampling to be near the grid scale
        centers = ilp(sam_pts, p, k, radius);
        // rescale back
        radius *= scale;
        for(auto &v: centers){
            v.x = v.x * scale + x_min;
            v.y = v.y * scale + y_min;
        }
    }
    printf("%.4lf\n", radius);
    for(auto v: centers){
        printf("%.4lf %.4lf\n", v.x, v.y);
    }
}