/*

Copyright (c) 2013  Ghassen Hamrouni

Abstract:

This module converts a bitmap image to a graph representation.

Author:

Ghassen Hamrouni <ghamrouni.iptech@gmail.com> 23-07-2013

Revision History:

*/

#pragma once

#ifndef ___IMAGE_GRAPH_H__
#define ___IMAGE_GRAPH_H__

#include <iostream>
#include <vector>

#include "Image.hpp"

#include <lemon/core.h>
#include <lemon/list_graph.h>
#include <lemon/smart_graph.h>
#include <lemon/concepts/maps.h>

#include <unordered_set>
#include <unordered_map>
#include <cmath>

struct Pixel {
	int x, y;
	int node;

	int getX() const {
		return x;
	}

	int getY() const {
		return y;
	}
};

class PixelHashFunction {
public:
	size_t operator ()(const Pixel &pix) const
	{
		int hashCode = 0;

		hashCode = hashCode * 31 + pix.x;
		hashCode = hashCode * 31 + pix.y;

		return hashCode;
	}
};

class PixelEqual {
public:
	bool operator ()(const Pixel &pix, const Pixel &pix2) const
	{
		return (pix.getX() == pix2.getX()) && (pix.getY() == pix2.getY());
	}
};

class ImageGraph {

protected:

	std::unordered_map<int, int> image2node;

public:

	ImageGraph() {
	}

	bool inRange(int x, int y, const Image& image) {
		return ((x >= 0) && (y >= 0) && (x < image.getWidth()) && (y < image.getHeight()));
	}

	Color getColor(int x, int y, const Image& image) {
		Color color;
		if ((x >= 0) && (y >= 0) && (x < image.getWidth()) && (y < image.getHeight())) {
			image.getPixel(x, y, color);
		}

		return color;
	}

	bool hasNode(const Image& image, int x, int y) {

		int imId = getNodeImageId(image, x, y);
		return (image2node.find(imId) != image2node.end());
	}

	int getNodeId(const Image& image, int x, int y) {
		if (hasNode(image, x, y)) {
			int imId = getNodeImageId(image, x, y);
			return image2node[imId];
		}

		return 0;
	}

	int getNodeImageId(const Image& image, int x, int y) {
		return y * image.getWidth() + x;
	}

	int count_nodes(const Image& image, float th) {
		int nb_nodes = 0;
		for (int x = 0; x < image.getWidth(); x++)
		{
			for (int y = 0; y < image.getHeight(); y++)
			{
				Color v = getColor(x, y, image);

				if (v.intensity() < th) {
					int imId = getNodeImageId(image, x, y);
					image2node[imId] = nb_nodes;
					nb_nodes += 1;
				}
			}
		}

		return nb_nodes;
	}

	// Convert an image to a graph using an 8-neighborhood system
	void build(const Image& image, R1::Graph* rgraph, float th)
	{
		int nb_nodes = count_nodes(image, th);
		rgraph->setSize(nb_nodes);
		int s = 50;

		for (int x = 0; x < image.getWidth(); x++)
		{
			for (int y = 0; y < image.getHeight(); y++)
			{
				int imId = getNodeImageId(image, x, y);
				
				//
				// 8-neighborhood system (Moore graph)
				// 
				Color v = getColor(x, y, image);

				if (v.intensity() < th) {

					auto node_1 = image2node[imId];

					if (inRange(x + 1, y, image))
					{
						Color w = getColor(x + 1, y, image);

						if (w.intensity() < th) {
							int w_id = getNodeImageId(image, x + 1, y);
							auto node_2 = image2node[w_id];

							if (!rgraph->is_connected(node_1, node_2)) {
								rgraph->connect(node_1, node_2, v.difference(w, s));
							}
						}
					}

					if (inRange(x, y + 1, image))
					{
						Color w = getColor(x, y + 1, image);

						if (w.intensity() < th) {
							int w_id = getNodeImageId(image, x, y + 1);
							auto node_2 = image2node[w_id];

							if (!rgraph->is_connected(node_1, node_2)) {
								rgraph->connect(node_1, node_2, v.difference(w, s));
							}
						}
					}

					if (inRange(x + 1, y + 1, image))
					{
						Color w = getColor(x + 1, y + 1, image);

						if (w.intensity() < th) {
							int w_id = getNodeImageId(image, x + 1, y + 1);
							auto node_2 = image2node[w_id];

							if (!rgraph->is_connected(node_1, node_2)) {
								rgraph->connect(node_1, node_2, v.difference(w, s));
							}
						}
					}

					if (inRange(x - 1, y + 1, image))
					{
						Color w = getColor(x - 1, y + 1, image);
						if (w.intensity() < th) {
							int w_id = getNodeImageId(image, x - 1, y + 1);
							auto node_2 = image2node[w_id];

							if (!rgraph->is_connected(node_1, node_2)) {
								rgraph->connect(node_1, node_2, v.difference(w, s));
							}
						}
					}

				}

			}
		}
	}

};

#endif // ___IMAGE_GRAPH_H__