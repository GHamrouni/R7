/*

Copyright (c) 2013  Ghassen Hamrouni

Abstract:

This module provides algorithms for weighted-graph clustering

Author:

Ghassen Hamrouni <ghamrouni.iptech@gmail.com> 11-04-2013

Revision History:

*/


#pragma once

#ifndef __NLP_CLUSTER_COLLECTION_H__
#define __NLP_CLUSTER_COLLECTION_H__

#include <memory>
#include <vector>

namespace R1
{
	// A node in a graph with
	// its associated cluster
	struct ClusteredNode
	{
		int id;
		int cluster_id;
	};

	struct ClusterCollection
	{
		int size;
		std::vector<ClusteredNode> clusters;
	};
}

#endif // __NLP_CLUSTER_COLLECTION_H__