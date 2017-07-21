/*

Copyright (c) 2013  Ghassen Hamrouni

Abstract:

This module provides a weighted undirected graph G (possibly with loops)

Author:

Ghassen Hamrouni <ghamrouni.iptech@gmail.com> 23-07-2013

Revision History:

*/

#ifndef R_1_PUBLIC_GRAPH_H_
#define R_1_PUBLIC_GRAPH_H_

#include <Eigen/Dense>
#include "MinCutClusterAnalyzer.hpp"
#include "ClusterCollection.hpp"

#include <random>
#include <fstream>
#include <string>
#include <set>
#include <unordered_set>

namespace R1 {

	//
	// Graph: A weighted undirected graph G (possibly with loops)
	//
	class Graph {

	public:

		virtual ~Graph() {}

		Graph() {}

		Graph(int NVertex) : n_vertex(NVertex), n_edges(0) {
			adjacency_matrix = Eigen::MatrixXf(NVertex, NVertex);
			adjacency_matrix.setZero();
		}

		void setSize(int NVertex) {
			n_vertex = (NVertex);
			n_edges  = (0);
			adjacency_matrix = Eigen::MatrixXf(NVertex, NVertex);
			adjacency_matrix.setZero();
		}

		void connect(int u, int v) {
			if (u >= n_vertex) {
				return;
			}

			if (v >= n_vertex) {
				return;
			}

			adjacency_matrix(u, v) = 1.0;
			adjacency_matrix(v, u) = 1.0;

			n_edges = n_edges + 1;
		}

		void connect(int u, int v, float weight) {
			if (u >= n_vertex) {
				return;
			}

			if (v >= n_vertex) {
				return;
			}

			adjacency_matrix(u, v) = weight;
			adjacency_matrix(v, u) = weight;

			n_edges = n_edges + 1;
		}

		float degree(int v) {
			float n = 0.0f;

			if (v >= n_vertex) {
				return 0;
			}

			for (int u = 0; u < n_vertex; u++) {
				n += adjacency_matrix(u, v);
			}

			return n;
		}

		bool is_connected(int u, int v) {
			if (u >= n_vertex) {
				return false;
			}

			if (v >= n_vertex) {
				return false;
			}

			return (adjacency_matrix(u, v) >= 1.0);
		}

		// A graph is said to be nontrivial if it contains at least
		// one edge.
		bool is_trivial() {
			return n_edges == 0;
		}

		void generateRandomSpanningTree()
		{
			std::random_device rd;
			std::mt19937 gen(rd());
			std::uniform_int_distribution<> vdist(0, n_vertex - 1);
			int current_vertex = vdist(gen);

			std::vector<int> P;
			std::vector<int> V;

			for (int i = 0; i < n_vertex; i++)
			{
				if (i != current_vertex)
					P.push_back(i);
			}

			V.push_back(current_vertex);

			while (!P.empty())
			{
				std::uniform_int_distribution<> distP(0, P.size() - 1);

				int j_index = distP(gen);
				int j = P[j_index];

				std::swap(P[j_index], P.back());
				P.pop_back();

				std::uniform_int_distribution<> distV(0, V.size() - 1);

				int u_index = distV(gen);
				int u = V[u_index];

				connect(u, j);
				
				V.push_back(j);
			}
		}

		void generateUniformRandomSpanningTree()
		{
			std::random_device rd;
			std::mt19937 gen(rd());
			std::uniform_int_distribution<> vdist(0, n_vertex - 1);

			int current_vertex = vdist(gen);
			std::vector<int> V;

			int n_edge = n_vertex - 1;

			for (int i = 0; i < n_vertex; i++)
			{
				V.push_back(0);
			}

			V[current_vertex] = 1;

			while (n_edge)
			{
				int j = vdist(gen);

				if (V[j] == 0)
				{
					V[j] = 1;
					connect(current_vertex, j);
					n_edge--;
				}

				current_vertex = j;
			}
		}

		void generateRandomGraph(int n_random_edges)
		{
			std::random_device rd;
			std::mt19937 gen(rd());
			std::uniform_int_distribution<> vdist(0, n_vertex - 1);
			
			while (n_random_edges > 0)
			{
				int i = vdist(gen);
				int j = vdist(gen);

				connect(i, j);

				n_random_edges--;
			}
		}

		ClusterCollection cluster(double alpha = 0.5)
		{
			return MinCutClusterAnalyzer::cluster(adjacency_matrix, alpha);
		}

		std::vector<ClusterCollection> hierarchicalCluster(double alpha = 0.5, double decay = 0.5)
		{
			std::vector<ClusterCollection> clusters;
			int prev_cluster_nb = 0;

			for (int i = 0; ; i++)
			{
				ClusterCollection cl = MinCutClusterAnalyzer::cluster(adjacency_matrix, alpha);

				if (cl.size != prev_cluster_nb)
				{
					prev_cluster_nb = cl.size;
					clusters.push_back(cl);
				}

				if (cl.size == 1) // One cluster
					break;

				if (alpha <= 0.0)
					break;

				alpha *= decay;
			}

			return clusters;
		}

	protected:
		Eigen::MatrixXf adjacency_matrix;

		int		n_vertex;
		int		n_edges;

	};
}
#endif