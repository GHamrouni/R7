/*

Copyright (c) 2013  Ghassen Hamrouni

Abstract:

This module provides algorithms for extracting topological minor

Author:

Ghassen Hamrouni <ghamrouni.iptech@gmail.com> 16-04-2013

Revision History:

*/


#pragma once

#ifndef __TOPOLOGICAL_MINOR_H__
#define __TOPOLOGICAL_MINOR_H__

#include <memory>

#include <lemon/core.h>
#include <lemon/list_graph.h>
#include <lemon/smart_graph.h>

namespace R1
{
	class TopologicalMinor
	{
	public:

		void extractSmallestMinor(lemon::ListGraph& graph, int maxIteration = 1000) {
			while (!extractMinor(graph) && (maxIteration > 0))
			{
				maxIteration--;
			}
		}

		bool extractMinor(lemon::ListGraph& graph) {

			bool invariant = true;

			for (lemon::ListGraph::NodeIt node(graph); node != lemon::INVALID; ++node)
			{
				int degree = 0;
				int n2 = 0;
				int n2prev = 0;

				int node_id = graph.id(node);

				for (lemon::ListGraph::IncEdgeIt e(graph, node); e != lemon::INVALID; ++e)
				{
					n2 = graph.id(graph.u(e));

					if (n2 == node_id)
						n2 = graph.id(graph.v(e));

					if (degree == 0) 
						n2prev = n2;

					degree++;
				}

				if (degree == 2) {

					auto node_2 = graph.nodeFromId(n2);
					bool valid = true;

					// Check if the removal of n2 creates a multi-graph
					for (lemon::ListGraph::IncEdgeIt e2(graph, node_2); e2 != lemon::INVALID; ++e2)
					{
						int incNodeId = graph.id(graph.u(e2));

						if (incNodeId == n2)
							incNodeId = graph.id(graph.v(e2));

						if (incNodeId == n2prev)
						{
							valid = false;
							break;
						}
					}

					if (valid)
					{
						graph.contract(node, node_2);
						invariant = false;
					}
				}
			}

			return invariant;
		}
	};
}

#endif // __TOPOLOGICAL_MINOR_H__