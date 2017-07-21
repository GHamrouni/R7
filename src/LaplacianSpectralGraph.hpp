/*

Copyright (c) 2013  Ghassen Hamrouni

Abstract:

This module provides some spectral graph algorithms

Author:

Ghassen Hamrouni <ghamrouni.iptech@gmail.com> 23-07-2013

Revision History:

*/

#ifndef R_1_PUBLIC_GRAPH_LAPLACIAN_SPECTRAL_H_
#define R_1_PUBLIC_GRAPH_LAPLACIAN_SPECTRAL_H_

#include <Eigen/Dense>
#include "SpectralGraph.hpp"

namespace R1 {
	class LaplacianSpectralGraph : public SpectralGraph {
	public:
		virtual ~LaplacianSpectralGraph() {}

		LaplacianSpectralGraph(int NVertex) : SpectralGraph(NVertex) {

		}

		//////////////////////////////////////////////////////////////////////////
		//								Laplacian
		//////////////////////////////////////////////////////////////////////////

		double laplacian(int u, int v) {

			double dv = degree(v);

			if (u == v && dv != 0) {
				return 1 - adjacency_matrix(v, v) / dv;
			}

			if (adjacency_matrix(u, v) >= 1.0) {

				double du = degree(u);

				return -adjacency_matrix(u, v) / sqrt(du * dv);
			}

			return 0.0;
		}

		Eigen::MatrixXd getLaplacianMatrix() {
			Eigen::MatrixXd laplacian_matrix(n_vertex, n_vertex);

			for (int i = 0; i < n_vertex; i++)	{
				for (int j = 0; j < n_vertex; j++)	{
					laplacian_matrix(i, j) = laplacian(i, j);
				}
			}

			return laplacian_matrix;
		}

		Eigen::VectorXd laplacianSpectrum() {

			Eigen::MatrixXd L = getLaplacianMatrix();

			Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eg1(L);
			return eg1.eigenvalues();
		}

		static bool isospectralLaplacian(LaplacianSpectralGraph g1, LaplacianSpectralGraph g2, double epsilon) {

			Eigen::VectorXd normalizedEg1 = g1.laplacianSpectrum().normalized();
			Eigen::VectorXd normalizedEg2 = g2.laplacianSpectrum().normalized();	

			if (normalizedEg1.cols() != normalizedEg2.cols()) {
				return false;
			}

			if (normalizedEg1.rows() != normalizedEg2.rows()) {
				return false;
			}

			for (int i = 0; i < normalizedEg1.rows(); i++)	{
				for (int j = 0; j < normalizedEg1.cols(); j++)	{
					if (abs(normalizedEg1(i, j) - normalizedEg2(i, j)) > epsilon) {
						return false;
					}
				}
			}

			return true;
		}

		static double distance(LaplacianSpectralGraph g1, LaplacianSpectralGraph g2) {
			Eigen::VectorXd normalizedEg1 = g1.laplacianSpectrum().normalized();
			Eigen::VectorXd normalizedEg2 = g2.laplacianSpectrum().normalized();	

			double dist = 0.0;

			if (normalizedEg1.cols() != normalizedEg2.cols()) {
				return 0.0;
			}

			if (normalizedEg1.rows() != normalizedEg2.rows()) {
				return 0.0;
			}

			for (int i = 0; i < normalizedEg1.rows(); i++)	{
				for (int j = 0; j < normalizedEg1.cols(); j++)	{
					dist += abs(normalizedEg1(i, j) - normalizedEg2(i, j));
				}
			}

			return dist;
		}
	};
}

#endif