/*

Copyright (c) 2013  G.H

Abstract:

This module provides algorithms for weighted-graph clustering

Author:

Ghassen Hamrouni <ghamrouni.iptech@gmail.com> 11-04-2013

Revision History:

*/


#pragma once

#ifndef __NLP_GRAPH_TOOLS_I_CLUSTERING_H__
#define __NLP_GRAPH_TOOLS_I_CLUSTERING_H__

#include <memory>
#include <vector>
#include <list>
#include <set>
#include <string>
#include <ctype.h>


#include <lemon/core.h>
#include <lemon/list_graph.h>
#include <lemon/smart_graph.h>
#include <lemon/concepts/maps.h>

#include "Eigen/Dense"

namespace R1
{
	class GraphConverter
	{
	public:

		typedef lemon::concepts::ReadMap<lemon::ListGraph::Edge, double> CapMap;
		typedef lemon::concepts::Graph Graph;

		// Transform an adjacency matrix to the corresponding lemon graph representation
		static void getGraph(lemon::ListGraph &graph, lemon::ListGraph::EdgeMap<double>& cmap, Eigen::MatrixXf adjacencyMatrix)
		{
			for (int i = 0; i < adjacencyMatrix.rows(); i++)
				graph.addNode();

			for (int i = 0; i < adjacencyMatrix.rows(); i++)
			{
				auto node_1 = lemon::ListGraph::nodeFromId(i);

				for (int j = 0; j < i; j++)
				{
					auto node_2 = lemon::ListGraph::nodeFromId(j);

					if (adjacencyMatrix(i, j) > 0.0)
					{
						auto e = graph.addEdge(node_1, node_2);
						cmap[e] = adjacencyMatrix(i, j);
					}
				}
			}
		}

		// Transform an adjacency matrix to the corresponding lemon graph representation
		static void getGraph(lemon::ListGraph &graph, Eigen::MatrixXf adjacencyMatrix)
		{
			for (int i = 0; i < adjacencyMatrix.rows(); i++)
				graph.addNode();

			for (int i = 0; i < adjacencyMatrix.rows(); i++)
			{
				auto node_1 = lemon::ListGraph::nodeFromId(i);

				for (int j = 0; j < i; j++)
				{
					auto node_2 = lemon::ListGraph::nodeFromId(j);

					if (adjacencyMatrix(i, j) > 0.0)
					{
						graph.addEdge(node_1, node_2);
					}
				}
			}
		}
	};
}

#endif // __NLP_GRAPH_TOOLS_I_CLUSTERING_H__