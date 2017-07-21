/*

Copyright (c) 2013  G.H

Abstract:

This module provides algorithms for weighted-graph clustering

Author:

Ghassen Hamrouni <ghamrouni.iptech@gmail.com> 11-04-2013

Revision History:

*/


#pragma once

#ifndef __NLP_RANDOM_WALK_CLUSTERING_H__
#define __NLP_RANDOM_WALK_CLUSTERING_H__

#include <memory>
#include <vector>
#include <list>
#include <set>
#include <string>
#include <ctype.h>


#include <lemon/core.h>
#include <lemon/list_graph.h>
#include <lemon/smart_graph.h>
#include <lemon/bfs.h>
#include <lemon/concepts/graph_components.h>
#include <lemon/connectivity.h>

#include "SpectralRank.hpp"
#include "GraphConverter.hpp"
#include "ClusterCollection.hpp"

namespace R1
{
	// Based on the paper:
	// On Clustering Using Random Walks David Harel and Yehuda Koren
	// Dept. of Computer Science and Applied Mathematics
	// The Weizmann Institute of Science, Rehovot, Israel

	class RWClusterAnalyzer
	{
	public:

		//
		// cluster : separate the graph-nodes in clusters
		//                               
		// Arguments:
		//      adjacencyMatrix		A symmetric-matrix representing the weighted
		//							undirected graph
		//      threshold			
		//		sIterations
		//		nRandomWalk
		//
		// Returns:
		//      Return the clustered nodes.
		//
		//
		static ClusterCollection cluster(
			Eigen::MatrixXf adjacencyMatrix, 
			double threshold,
			int sIterations, 
			int nRandomWalk)
		{
			// NS: Separation by neighborhood similarity
			auto NS = calculateSeparation(adjacencyMatrix, nRandomWalk, sIterations);

			//threshold = calculateThreshold(adjacencyMatrix, NS); // TODO: test this
			auto TH = thresholdMatrix(NS, threshold);

			ClusterCollection clusterCollection;

			lemon::ListGraph graph;
			GraphConverter::getGraph(graph, TH);

			lemon::ListGraph::NodeMap<int> nmap(graph);
			int nc = lemon::connectedComponents(graph, nmap);

			clusterCollection.size = nc;

			// Get the connected components from nmap
			for (int i = 0; i < adjacencyMatrix.rows(); i++)
			{
				auto node = lemon::ListGraph::nodeFromId(i);

				ClusteredNode cnode;
				cnode.id = i;
				cnode.cluster_id = nmap[node];

				clusterCollection.clusters.push_back(cnode);
			}

			return clusterCollection;
		}

		// Given a weighted graph returns the associated transition matrix
		static Eigen::MatrixXf getTransitionMatrix(Eigen::MatrixXf adjacencyMatrix)
		{
			Eigen::MatrixXf Mg = SpectralRanking::SetStochastic(adjacencyMatrix);

			return Mg;
		}

		// Calculate the transition matrix after K iterations
		static Eigen::MatrixXf performKVisit(Eigen::MatrixXf transitionMatrix, int k)
		{
			Eigen::MatrixXf T = transitionMatrix;

			for (int i = 0; i < k; ++i)
			{
				T = T * transitionMatrix;
			}

			return T;
		}

		// Sum all the transition matrices from 0 to K iteration
		static Eigen::MatrixXf performKVisitSum(Eigen::MatrixXf transitionMatrix, int k)
		{
			Eigen::MatrixXf T = transitionMatrix;
			Eigen::MatrixXf S = T;
			double N = k;

			for (int i = 0; i < k; ++i)
			{
				T = T * transitionMatrix;
				S += T;
			}

			return S;
		} 

		// Get the vector whose j-th component is the probability that a random walk
		// originating at i will visit node j in its k-th step
		static Eigen::VectorXf getTransitionVector(int sourceNode, Eigen::MatrixXf transitionMatrix, int k)
		{
			Eigen::MatrixXf T = transitionMatrix;

			for (int i = 0; i < k; ++i)
			{
				T = T * transitionMatrix;
			}

			return T.row(sourceNode);
		}

		// Sum all the transition vectors from 0 to K iteration
		static Eigen::VectorXf getTransitionVectorSum(int sourceNode, Eigen::MatrixXf transitionMatrix, int k)
		{
			Eigen::MatrixXf T = transitionMatrix;
			Eigen::VectorXf tr_i = T.row(sourceNode);

			for (int i = 0; i < k; ++i)
			{
				T = T * transitionMatrix;
				tr_i += T.row(sourceNode);
			}

			return tr_i;
		}

		// A simple similarity measure between two vectors
		// large values = similar vectors
		static float similarityMeasure(Eigen::VectorXf x1, Eigen::VectorXf x2, int k)
		{
			float lp = 0.0;

			if (x1.maxCoeff() != 0.0 || x2.maxCoeff() != 0.0)
			{
				if (x1.maxCoeff() == 0.0 || x2.maxCoeff() == 0.0)
				{
					return 0.0;
				}
			}

			for (int i = 0; i < x1.size(); ++i)
			{
				lp += fabs(x1(i) - x2(i));
			}

			float fk = fabs(expf(2 * k - lp) - 1.0f);

			return fk;
		}

		// NS: Separation by neighborhood similarity
		static Eigen::MatrixXf calculateSeparation(Eigen::MatrixXf adjacencyMatrix, int k)
		{
			Eigen::MatrixXf transitionMatrix = getTransitionMatrix(adjacencyMatrix);
			Eigen::MatrixXf sumTransitions   = performKVisitSum(transitionMatrix, k);

			Eigen::MatrixXf separationMatrix(adjacencyMatrix.rows(), adjacencyMatrix.cols());
			separationMatrix.fill(0.0);

			for (int i = 0; i < adjacencyMatrix.rows(); i++)
			{
				for (int j = 0; j < adjacencyMatrix.cols(); j++)
				{
					if (i != j)
					{
						Eigen::VectorXf x1 = sumTransitions.row(i);
						Eigen::VectorXf x2 = sumTransitions.row(j);

						float ws = similarityMeasure(x1, x2, k);

						separationMatrix(i, j) = ws;
						separationMatrix(j, i) = ws;
					}
				}
			}

			return separationMatrix;
		}

		// NS: Separation by neighborhood similarity
		static Eigen::MatrixXf calculateSeparation(Eigen::MatrixXf adjacencyMatrix, int k, int iterations)
		{
			Eigen::MatrixXf S = adjacencyMatrix;

			for (int i = 0; i < iterations; i++)
			{
				S = calculateSeparation(S, k);
			}

			return S;
		}

#define HISTO_SIZE 10

		// Try to calculate a threshold
		static double calculateThreshold(Eigen::MatrixXf& adjacencyMatrix, Eigen::MatrixXf& transitionMatrix)
		{
			int H[HISTO_SIZE];

			double M = transitionMatrix.maxCoeff();

			for (int i = 0; i < HISTO_SIZE; ++i)
			{
				H[i] = 0;
			}

			if (M == 0.0)
				return 0.0;

			for (int i = 0; i < adjacencyMatrix.rows(); i++)
			{
				for (int j = 0; j < adjacencyMatrix.cols(); j++)
				{
					if (adjacencyMatrix(i, j) > 0)
					{
						int u = (int) ((HISTO_SIZE - 1) * (transitionMatrix(i, j) / M));

						H[u] += 1;
					}
				}
			}

			double minHindex = -1;
			double maxHindex = H[0];

			for (int i = 0; i < HISTO_SIZE; ++i)
			{
				if (H[i] > maxHindex)
					maxHindex = i;

				if (H[i] != 0 && (minHindex == -1 || H[i] < minHindex))
					minHindex = i;
			}

			if (minHindex == -1)
				minHindex = 0;

			double index = (minHindex + maxHindex) / 2.0;

			double TH = index * M / double(HISTO_SIZE);
			return TH;
		}

		// Set values  < threshold to 0
		// and values >= threshold to 1
		static Eigen::MatrixXf thresholdMatrix(Eigen::MatrixXf m, double v)
		{
			Eigen::MatrixXf thM(m.rows(), m.cols());
			thM.fill(0.0);

			for (int i = 0; i < m.rows(); i++)
			{
				for (int j = 0; j < m.cols(); j++)
				{
					if (m(i, j) >= v)
					{
						thM(i, j) = 1.0;
					}
				}
			}

			return thM;
		}
	};
}

#endif // __NLP_RANDOM_WALK_CLUSTERING_H__