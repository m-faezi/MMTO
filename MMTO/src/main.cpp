#include "pybind11/pybind11.h"
#include <pybind11/stl.h>
#include "higra/graph.hpp"
#include "higra/algo/tree.hpp"
#include "higra/attribute/tree_attribute.hpp"
#include <xtensor/xnoalias.hpp>
#include <vector>
#include <stack>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
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

            // Create CSV content
            std::stringstream csv_content;
            csv_content << "tree_i_id,object_i_id,tree_j_id,object_j_id,flux_i,flux_j,cosine_similarity,distance\n";

            int match_count = 0;

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

                        j < (index_t) ntrees && tj != ti;
                        tj++, j++, mj++, xj++, yj++, dj++, aj++, area_j++, moment_j++
                    ) {

                        (**tj).compute_children();
                        auto ses_ij_n = ses(i, j)(n);

                        double dist = sqrt(
                            ((**xi)[n] - (**xj)[ses_ij_n])*((**xi)[n] - (**xj)[ses_ij_n]) +
                            ((**yi)[n] - (**yj)[ses_ij_n])*((**yi)[n] - (**yj)[ses_ij_n])
                        );

                        if (dist < 3) {

                            double X[] = {
                                (**area_i)[n],
                                (**moment_i)[n],
                                (**mi)[n],
                                (**mi)[n]/(**area_i)[n]
                            };
                            double Y[] = {
                                (**area_j)[ses_ij_n],
                                (**moment_j)[ses_ij_n],
                                (**mj)[ses_ij_n],
                                (**mj)[ses_ij_n]/(**area_j)[ses_ij_n]
                            };

                            double cos_sim = cosine_similarity(X, Y, 4);

                            if (cos_sim > 0.3) {
                                // Add row to CSV
                                csv_content << i << "," << n << ","
                                           << j << "," << ses_ij_n << ","
                                           << (**mi)[n] << "," << (**mj)[ses_ij_n] << ","
                                           << cos_sim << "," << dist << "\n";
                                match_count++;
                            }
                        }
                    }
                }
            }

            // Write CSV to file
            std::ofstream csv_file("detection_colors.csv");
            if (csv_file.is_open()) {
                csv_file << csv_content.str();
                csv_file.close();
                std::cout << "CSV file 'object_matches.csv' created with " << match_count << " matches." << std::endl;
            } else {
                std::cerr << "Error: Could not create CSV file." << std::endl;
            }

            // Return the CSV content as string (optional)
            return csv_content.str();
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

