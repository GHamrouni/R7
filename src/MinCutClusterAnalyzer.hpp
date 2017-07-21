/*

Copyright (c) 2013  G.H

Abstract:

This module provides algorithms for weighted-graph clustering

Author:

Ghassen Hamrouni <ghamrouni.iptech@gmail.com> 16-04-2013

Revision History:

*/


#pragma once

#ifndef __NLP_MIN_CUT_CLUSTERING_H__
#define __NLP_MIN_CUT_CLUSTERING_H__

#include <memory>
#include <vector>
#include <list>
#include <set>
#include <string>
#include <ctype.h>
#include <unordered_map>

#include <lemon/core.h>
#include <lemon/list_graph.h>
#include <lemon/smart_graph.h>
#include <lemon/bfs.h>
#include <lemon/concepts/graph_components.h>
#include <lemon/concepts/maps.h>
#include <lemon/connectivity.h>
#include <lemon/gomory_hu.h>


#include "GraphConverter.hpp"
#include "ClusterCollection.hpp"

namespace R1
{
	// Based on Graph Clustering and Minimum Cut Trees
	// Gary William Flake, Robert E. Tarjan, and Kostas Tsioutsiouliklis
	class MinCutClusterAnalyzer
	{
	public:

		//
		// cluster : separate the graph-nodes in clusters
		//                               
		// Arguments:
		//      adjacencyMatrix		A symmetric-matrix representing the weighted
		//							undirected graph
		//		alpha				Controls the number of clusters
		//
		// Returns:
		//      Return the clustered nodes.
		//
		//
		static ClusterCollection cluster(lemon::ListGraph& graph, lemon::ListGraph::EdgeMap<double>& capacity, double alpha = 0.5)
		{
			typedef lemon::ListGraph::EdgeMap<double> CapMap;
			typedef lemon::GomoryHu<lemon::ListGraph, CapMap> GomoryHuT;

			std::unordered_map<int, int>	nodeClusterMap; // Node Id --> Cluster Id

			ClusterCollection clusterCollection;

			lemon::ListGraph::NodeMap<bool> cutMap(graph);

			// Add an artificial sink t
			auto t = graph.addNode();

			for (lemon::ListGraph::NodeIt n(graph); n != lemon::INVALID; ++n)
			{
				auto e = graph.addEdge(n, t);
				capacity[e] = alpha;
			}

			lemon::GomoryHu<lemon::ListGraph, CapMap> ght(graph, capacity);

			ght.run(); // Build the Gomory-Hu Tree

			int nc = 0; // The number of different clusters

			// Remove the artificial sink t from the min-cut tree, and get 
			// the resulting connected components form the clusters of G
			for (lemon::ListGraph::NodeIt node(graph); node != lemon::INVALID; ++node)
			{
				if (node == t) 
					continue;

				std::list<int> currentBranch;
				int id = lemon::ListGraph::id(node);

				// Verify if the node is already clustered
				if (nodeClusterMap.find(id) == nodeClusterMap.end())
				{
					currentBranch.push_back(id);

					auto p = ght.predNode(node);
					int  prev_node_id = lemon::ListGraph::id(p);

					while (p != t && prev_node_id != -1 && 
						(nodeClusterMap.find(prev_node_id) == nodeClusterMap.end()))
					{
						currentBranch.push_back(prev_node_id);
						p = ght.predNode(p);

						prev_node_id = lemon::ListGraph::id(p);
					}

					int clusterId = nc;

					if (nodeClusterMap.find(prev_node_id) == nodeClusterMap.end())
					{
						nc++;
					} else {
						clusterId = nodeClusterMap[prev_node_id]; // The branch is a subbranch
					}

					// Assign a cluster to each node in the current branch
					for (auto it = currentBranch.begin(); it != currentBranch.end(); ++it)
					{							
						ClusteredNode cnode;
						cnode.id = *it;
						cnode.cluster_id = clusterId;

						clusterCollection.clusters.push_back(cnode);

						nodeClusterMap[*it] = clusterId;
					}
				}
			}

			clusterCollection.size = nc;

			// Remove the sink node
			graph.erase(t);

			return clusterCollection;
		}

		static std::vector<ClusterCollection> hierarchicalCluster(
			lemon::ListGraph& graph, 
			lemon::ListGraph::EdgeMap<double>& capacity, 
			double alpha = 0.5, 
			double decay = 0.5)
		{
			std::vector<ClusterCollection> clusters;
			int prev_cluster_nb = 0;

			for (int i = 0; ; i++)
			{
				ClusterCollection cluster = MinCutClusterAnalyzer::cluster(graph, capacity, alpha);

				if (cluster.size != prev_cluster_nb)
				{
					prev_cluster_nb = cluster.size;
					clusters.push_back(cluster);
				}

				if (cluster.size == 1) // One cluster
					break;

				if (alpha <= 0.0)
					break;

				alpha *= decay;
			}

			return clusters;
		}

		//
		// cluster : separate the graph-nodes in clusters
		//                               
		// Arguments:
		//      adjacencyMatrix		A symmetric-matrix representing the weighted
		//							undirected graph
		//		alpha				Controls the number of clusters
		//
		// Returns:
		//      Return the clustered nodes.
		//
		//
		static ClusterCollection cluster(const Eigen::MatrixXf& adjacencyMatrix, double alpha = 0.5)
		{
			typedef lemon::ListGraph::EdgeMap<double> CapMap;

			lemon::ListGraph graph;
			CapMap capacity(graph);

			GraphConverter::getGraph(graph, capacity, adjacencyMatrix);


			return cluster(graph, capacity, alpha);
		}
	};
}

#endif // __NLP_MIN_CUT_CLUSTERING_H__