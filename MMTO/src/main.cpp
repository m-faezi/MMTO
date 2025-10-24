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
            const std::vector<hg::tree> & trees,
            const std::vector<xt::pyarray<double>> & xs,
            const std::vector<xt::pyarray<double>> & ys,
            const std::vector<xt::pyarray<double>> & fluxes,
            const std::vector<xt::pyarray<double>> & gammas,
            const std::vector<xt::pyarray<double>> & as,
            const std::vector<xt::pyarray<double>> & volumes,
            const std::vector<xt::pyarray<double>> & ids,
            const std::vector<std::string> & tree_ids
        ) {

            auto ntrees = trees.size();

            if (ntrees == 0) {
                return std::string(""); // Return empty string if no trees
            }

            auto nleaves = num_leaves(trees[0]);

            vector<vector<index_t>> adj_lists;
            vector<array_1d<index_t>> areas;
            array_2d<array_1d<index_t>> ses = xt::empty<array_1d<index_t>>({ntrees, ntrees});

            for (index_t i = 0; i < ntrees; i++) {
                areas.push_back(attribute_area(trees[i]));
                for (index_t j = 0; j < ntrees; j++) {
                    if (j != i) {
                        ses(i, j) = attribute_smallest_enclosing_shape(trees[i], trees[j]);
                    }
                }
            }

            vector<array_1d<index_t>> node_maps;
            adj_lists.resize(nleaves);

            for (index_t i = 0; i < ntrees; i++) {
                node_maps.emplace_back(array_1d<index_t>::from_shape({num_vertices(trees[i])}));
                xt::noalias(xt::view(node_maps[i], xt::range(0, nleaves))) = xt::arange<index_t>(nleaves);

                for (index_t n: leaves_to_root_iterator(trees[i], leaves_it::exclude, root_it::exclude)) {
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

            for (index_t i = 0; i < ntrees; i++) {
                node_maps[i](root(trees[i])) = rootn;
            }

            std::stringstream csv_content;
            csv_content << "tree_i_id,object_i_id,tree_j_id,object_j_id,flux_i,flux_j,cosine_similarity,distance\n";

            int match_count = 0;

            for (index_t i = 0; i < ntrees; i++) {
                for (index_t n: leaves_to_root_iterator(trees[i], leaves_it::include, root_it::exclude)) {
                    auto represent_n = node_maps[i](n);

                    // APPLY FILTER ONLY HERE - check both nodes have valid segment IDs
                    if (ids[i][n] <= 0) {
                        continue;
                    }

                    adj_lists[node_maps[i](parent(n, trees[i]))].push_back(represent_n);

                    for (index_t j = 0; j < ntrees; j++) {
                        if (j <= i) continue;

                        trees[j].compute_children();
                        auto ses_ij_n = ses(i, j)(n);

                        // APPLY FILTER ONLY HERE - check both nodes have valid segment IDs
                        if (ids[j][ses_ij_n] <= 0) {
                            continue;
                        }

                        double xi_val = xs[i][n];
                        double xj_val = xs[j][ses_ij_n];
                        double yi_val = ys[i][n];
                        double yj_val = ys[j][ses_ij_n];

                        double dist = sqrt(
                            (xi_val - xj_val) * (xi_val - xj_val) +
                            (yi_val - yj_val) * (yi_val - yj_val)
                        );

                        if (dist < 3) {

                            double area_i_val = as[i][n];
                            double volume_i_val = volumes[i][n];
                            double flux_i_val = fluxes[i][n];

                            double area_j_val = as[j][ses_ij_n];
                            double volume_j_val = volumes[j][ses_ij_n];
                            double flux_j_val = fluxes[j][ses_ij_n];

                            double X[] = {
                                area_i_val,
                                volume_i_val,
                                flux_i_val,
                                flux_i_val / area_i_val
                            };
                            double Y[] = {
                                area_j_val,
                                volume_j_val,
                                flux_j_val,
                                flux_j_val / area_j_val
                            };

                            double cos_sim = cosine_similarity(X, Y, 4);

                            double n_id = ids[i][n];
                            double m_id = ids[j][ses_ij_n];

                            if (cos_sim > 0.93) {
                                csv_content << tree_ids[i] << "," << n_id << ","
                                           << tree_ids[j] << "," << m_id << ","
                                           << flux_i_val << "," << flux_j_val << ","
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

