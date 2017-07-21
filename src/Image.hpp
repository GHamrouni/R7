/*

Copyright (c) 2013  Ghassen Hamrouni

Abstract:

This module provides IO functions for PNG images.

Author:

Ghassen Hamrouni <ghamrouni.iptech@gmail.com> 23-07-2013

Revision History:

*/

#pragma once

#ifndef ___IMAGE_PNG_H__
#define ___IMAGE_PNG_H__

#include <iostream>
#include <vector>
#include "lodepng.h"

#include "Color.hpp"

class Image {

	std::vector<unsigned char> buffer;
	unsigned int width, height;

public:

	int save(const char* filename)
	{
		unsigned error = lodepng::encode(filename, buffer, width, height);
		return error;
	}

	int load(const char* filename) {
		unsigned error = lodepng::decode(buffer, width, height, filename);	
		return error;
	}

	void init(unsigned int width, unsigned int height) {

		this->width = width;
		this->height = height;

		buffer.clear();
		
		for (unsigned int i = 0; i < width * height * 4; i++)
			buffer.push_back(255);
	}

	void setPixel(int x, int y, Color& color) {
		buffer[(y * width + x) * 4] = color.R;
		buffer[(y * width + x) * 4 + 1] = color.G;
		buffer[(y * width + x) * 4 + 2] = color.B;
		buffer[(y * width + x) * 4 + 3] = color.A;
	}

	void getPixel(int x, int y, Color& color) const {
		color.R = buffer[(y * width + x) * 4];
		color.G = buffer[(y * width + x) * 4 + 1];
		color.B = buffer[(y * width + x) * 4 + 2];
		color.A = buffer[(y * width + x) * 4 + 3];
	}

	int pixelCount() const {
		return width * height;
	}

	int getWidth() const {
		return width;
	}

	int getHeight() const {
		return height;
	}
};

#endif // ___IMAGE_PNG_H__