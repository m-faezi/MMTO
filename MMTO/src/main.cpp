#include "pybind11/pybind11.h"
#include <pybind11/stl.h>
#include "higra/graph.hpp"
#include "higra/algo/tree.hpp"   //to use is_leaf(node, tree)
#include "higra/attribute/tree_attribute.hpp"
#include <xtensor/xnoalias.hpp>
#include <vector>
#include <stack>
#include <iostream>
#define FORCE_IMPORT_ARRAY
#include "xtensor-python/pyarray.hpp"


namespace hg {

    namespace tree_fusion_internal {

        using namespace std;
        using index_t = int64_t;
        const index_t invalid_index = -1;

        float distance(int x1, int y1, int x2, int y2)
        {
            return sqrt(pow(x2 - x1, 2) + pow(y2 - y1, 2) * 1.0);
        }

        double correlation_coefficient(double X[], double Y[], double n) {

            double sum_X = 0, sum_Y = 0, sum_XY = 0;
            double squareSum_X = 0, squareSum_Y = 0;

            for (int i = 0; i < n; i++){

                sum_X = sum_X + X[i];
                sum_Y = sum_Y + Y[i];
                sum_XY = sum_XY + X[i] * Y[i];
                squareSum_X = squareSum_X + X[i] * X[i];
                squareSum_Y = squareSum_Y + Y[i] * Y[i];
            }

            double corr = (double)(n * sum_XY - sum_X * sum_Y)
                / sqrt((n * squareSum_X - sum_X * sum_X)
                * (n * squareSum_Y - sum_Y * sum_Y));

            return corr;
        }

        // Pearson
        double mean(int n, double arr[]) {
            double sum = 0;
            for (int i = 0; i < n; i++)
                sum += arr[i];

            return sum / n;
        }

        double stdDev(int n, double arr[], double mean) {
            double sum = 0;
            for (int i = 0; i < n; i++)
                sum += pow(abs(arr[i] - mean), 2);

            return sqrt(sum / n);
        }

        double pearson(double X[], double Y[], double n) {
            double xMean = mean(n, X);
            double yMean = mean(n, Y);
            double xStdDev = stdDev(n, X, xMean);
            double yStdDev = stdDev(n, Y, yMean);

            double sum = 0;
            for (int i = 0; i < n; i++)
                sum += (X[i] - xMean) * (Y[i] - yMean);

            return sum / (n * xStdDev * yStdDev);
        }

        double cosine_similarity(double A[], double B[], double n)
        {
            double dot = 0.0, denom_a = 0.0, denom_b = 0.0 ;
             for(int i = 0; i < n; i++) {
                dot += A[i] * B[i] ;
                denom_a += A[i] * A[i] ;
                denom_b += B[i] * B[i] ;
            }
            return dot / (sqrt(denom_a) * sqrt(denom_b)) ;
        }

        auto tree_map(
            const hg::tree & tree_1,
            const hg::tree & tree_2,
            const hg::tree & tree_3,
            const xt::pyarray<double> & mu_1,
            const xt::pyarray<double> & mu_2,
            const xt::pyarray<double> & mu_3,
            const xt::pyarray<double> & x_1,
            const xt::pyarray<double> & x_2,
            const xt::pyarray<double> & x_3,
            const xt::pyarray<double> & y_1,
            const xt::pyarray<double> & y_2,
            const xt::pyarray<double> & y_3,
            const xt::pyarray<double> & d_1,
            const xt::pyarray<double> & d_2,
            const xt::pyarray<double> & d_3,
            const xt::pyarray<double> & a_1,
            const xt::pyarray<double> & a_2,
            const xt::pyarray<double> & a_3,
            const xt::pyarray<double> & norm_area_1,
            const xt::pyarray<double> & norm_area_2,
            const xt::pyarray<double> & norm_area_3,
            const xt::pyarray<double> & moment_1,
            const xt::pyarray<double> & moment_2,
            const xt::pyarray<double> & moment_3
        ) {
            std::vector<const hg::tree *> trees{
                &tree_1,
                &tree_2,
                &tree_3
            };
            auto first = trees.begin(), last = trees.end();

            std::vector<const xt::pyarray<double> *> mus{
                &mu_1,
                &mu_2,
                &mu_3
            };
            auto first_mu = mus.begin(), last_mu = mus.end();

            std::vector<const xt::pyarray<double> *> xs{
                &x_1,
                &x_2,
                &x_3
            };
            auto first_x = xs.begin(), last_x = xs.end();

            std::vector<const xt::pyarray<double> *> ys{
                &y_1,
                &y_2,
                &y_3
            };
            auto first_y = ys.begin(), last_y = ys.end();

            std::vector<const xt::pyarray<double> *> ds{
                &d_1,
                &d_2,
                &d_3
            };
            auto first_d = ds.begin(), last_d = ds.end();

            std::vector<const xt::pyarray<double> *> as{
                &a_1,
                &a_2,
                &a_3
            };
            auto first_a = as.begin(), last_a = as.end();

            std::vector<const xt::pyarray<double> *> norm_areas{
                &norm_area_1,
                &norm_area_2,
                &norm_area_3
            };
            auto first_norm_area = norm_areas.begin(), last_norm_area = norm_areas.end();

            std::vector<const xt::pyarray<double> *> moments{
                &moment_1,
                &moment_2,
                &moment_3
            };
            auto first_moment = moments.begin(), last_moment = moments.end();

            index_t i, j;
            auto ti = first, tj = first;
            auto mi = first_mu, mj = first_mu;
            auto xi = first_x, xj = first_x;
            auto yi = first_y, yj = first_y;
            auto di = first_d, dj = first_d;
            auto ai = first_a, aj = first_a;
            auto area_i = first_norm_area, area_j = first_norm_area;
            auto moment_i = first_moment, moment_j = first_moment;
            auto ntrees = last - first;
            auto nleaves = num_leaves(**first);

            vector<vector<index_t>> adj_lists;
            vector<array_1d<index_t>> areas;
            array_2d<array_1d<index_t>> ses = xt::empty<array_1d<index_t>>({ntrees, ntrees});

            for (ti = first, i = 0; ti != last; ti++, i++) {
                areas.push_back(attribute_area(**ti));
                for (tj = first, j = 0; tj != last; tj++, j++) {
                    if (j != i) {
                        ses(i, j) = attribute_smallest_enclosing_shape(**ti, **tj);
                    }
                }
            }

            vector<array_1d<index_t>> node_maps;
            adj_lists.resize(nleaves);

            for (ti = first, i = 0; ti != last; ti++, i++) {
                node_maps.emplace_back(array_1d<index_t>::from_shape({num_vertices(**ti)}));
                xt::noalias(xt::view(node_maps[i], xt::range(0, nleaves))) = xt::arange<index_t>(nleaves);

                for (index_t n: leaves_to_root_iterator(**ti, leaves_it::exclude, root_it::exclude)) {
                    bool keep = true;

                    for (index_t j = 0; j < i && keep; j++) {
                        auto ses_ij_n = ses(i, j)(n);
                        if (areas[j](ses_ij_n) == areas[i](n)) {
                            keep = false;
                            node_maps[i](n) = node_maps[j](ses_ij_n);
                        }
                    }
                    if (keep) {
                        node_maps[i](n) = adj_lists.size();
                        adj_lists.emplace_back();
                    }
                }
            }

            auto rootn = adj_lists.size();
            adj_lists.emplace_back();

            for (ti = first, i = 0; ti != last; ti++, i++) {
                node_maps[i](root(**ti)) = rootn;
            }

            for (
                ti = first, mi = first_mu, xi = first_x, yi = first_y, di = first_d, ai = first_a,
                area_i = first_norm_area, moment_i = first_moment, i = 0;
                ti != last;
                ti++, i++, mi++, xi++, yi++, di++, ai++, area_i++, moment_i++
            ) {

                for (index_t n: leaves_to_root_iterator(**ti, leaves_it::include, root_it::exclude)) {
                    auto represent_n = node_maps[i](n);
                    adj_lists[node_maps[i](parent(n, **ti))].push_back(represent_n);

                    for (
                        tj = first, mj = first_mu, xj = first_x, yj = first_y, dj = first_d, aj = first_a,
                        area_j = first_norm_area, moment_j = first_moment, j = 0;
                        // j < (index_t) ntrees;
                        j < (index_t) ntrees && tj != ti;
                        tj++, j++, mj++, xj++, yj++, dj++, aj++, area_j++, moment_j++
                    ) {
                        (**tj).compute_children();
                        auto ses_ij_n = ses(i, j)(n);
                        if (
                            (i < j) && !((**xi)[n] < 0)
                        ){

                            if (
                                !((**xj)[ses_ij_n] < 0) &&
                                    (
                                        sqrt(
                                            ((**xi)[n] - (**xj)[ses_ij_n])*((**xi)[n] - (**xj)[ses_ij_n]) +
                                            ((**yi)[n] - (**yj)[ses_ij_n])*((**yi)[n] - (**yj)[ses_ij_n])
                                            //) < min(sqrt(areas[i](n))/3.14, sqrt(areas[j](ses_ij_n))/3.14)
                                        ) < 10
                                    )
                            ){

                                double X[] = {
                                    //areas[i](n),
                                    (**area_i)[n],
                                    (**moment_i)[n],
                                    (**mi)[n],
                                    (**mi)[n]/(**area_i)[n]
                                    //(**mi)[n]/areas[i](n)
                                };
                                double Y[] = {
                                    //areas[j](ses_ij_n),
                                    (**area_j)[ses_ij_n],
                                    (**moment_j)[ses_ij_n],
                                    (**mj)[ses_ij_n],
                                    (**mj)[ses_ij_n]/(**area_j)[ses_ij_n]
                                    //(**mj)[ses_ij_n]/areas[j](ses_ij_n)
                                };

                                if (cosine_similarity(X, Y, 4) > 0.5){
                                    //adj_lists[node_maps[j](ses_ij_n)].push_back(node_maps[i](n));
                                    adj_lists[node_maps[j](parent(ses_ij_n, **tj))].push_back(node_maps[i](n));
                                }
                            }


                            for (auto c_c: children_iterator(ses_ij_n, **tj)) {
                                if (
                                    !((**xj)[c_c] < 0) &&
                                    (
                                        sqrt(
                                            ((**xi)[n] - (**xj)[c_c])*((**xi)[n] - (**xj)[c_c]) +
                                            ((**yi)[n] - (**yj)[c_c])*((**yi)[n] - (**yj)[c_c])
                                        )
                                        // < min(sqrt(areas[i](n))/3.14, sqrt(areas[j](c_c))/3.14)
                                        < 10
                                    )
                                    // && ((**di)[n] == (**dj)[c_c])
                                ){

                                    double X[] = {
                                        //areas[i](n),
                                        (**area_i)[n],
                                        (**moment_i)[n],
                                        (**mi)[n],
                                        (**mi)[n]/(**area_i)[n]
                                        //(**mi)[n]/areas[i](n)
                                    };
                                    double Y[] = {
                                        //areas[j](c_c),
                                        (**area_j)[c_c],
                                        (**moment_j)[c_c],
                                        (**mj)[c_c],
                                        (**mj)[c_c]/(**area_j)[c_c]
                                        //(**mj)[c_c]/areas[j](c_c)
                                    };
                                    if (cosine_similarity(X, Y, 4) > 0.5){
                                        //adj_lists[node_maps[j](ses_ij_n)].push_back(node_maps[i](n));
                                        adj_lists[node_maps[j](parent(c_c, **tj))].push_back(node_maps[i](n));
                                    }

                                    else{
                                        //if (areas[j](c_c) <= areas[i](n)){
                                        if ((**area_j)[c_c] <= (**area_i)[n]){
                                            adj_lists[node_maps[i](n)].push_back(node_maps[j](c_c));
                                        }
                                        else{
                                            adj_lists[node_maps[j](c_c)].push_back(node_maps[i](n));
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }

            auto nnodes = adj_lists.size();
            array_1d<index_t> sorted_nodes = xt::empty<index_t>({nnodes});
            array_1d<char> marks = xt::zeros<char>({nnodes});
            stack<index_t> s;

            index_t count = 0;
            s.push(rootn);
            while (!s.empty()) {
                auto n = s.top();
                if (marks(n) > 0) {
                    s.pop();
                    if (marks(n) == 1) {
                        sorted_nodes(count++) = n;
                        marks(n) = 2;
                    }
                } else {
                    marks(n) = 1;
                    for (auto o: adj_lists[n]) {
                        if (marks(o) != 2) {
                            s.push(o);
                        }
                    }
                }
            }

            xt::pyarray<double> depth = xt::zeros<index_t>({nnodes});

            for (index_t i = nnodes - 1; i >= 0; i--) {
                index_t n = sorted_nodes[i];
                for (auto o: adj_lists[n]) {
                    depth(o) = (std::max)(depth(o), depth(n) + 1);
                }
            }

            return xt::eval(xt::view(depth, xt::range(0, nleaves)));

        }

        PYBIND11_MODULE(mmto, m)
        {
            xt::import_numpy();

            m.doc() = R"pbdoc(
                Multi-dimensional faint object detection
            )pbdoc";

            m.def("tree_map", tree_map, "Tree map function.");
        }
    }
}


