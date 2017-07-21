/*

Copyright (c) 2013  Ghassen Hamrouni

Abstract:

	This module provides a representation of an RGBA color

Author:

Ghassen Hamrouni <ghamrouni.iptech@gmail.com> 23-07-2013

Revision History:

*/

#pragma once

#ifndef ___RGBA_COLOR_IMAGING_H__
#define ___RGBA_COLOR_IMAGING_H__

static float epsilonXYZ = 0.008856f;
static float kXYZ = 903.3f;

#include <math.h>

struct Color {
	
	unsigned char R;
	unsigned char G;
	unsigned char B;
	unsigned char A;

	Color(	
		unsigned char r,
		unsigned char g,
		unsigned char b,
		unsigned char a) 
	{
		R = r;
		G = g;
		B = b;
		A = a;
	}

	Color() {
		R = 0;
		G = 0;
		B = 0;
		A = 0;
	}

	static float linear(float V) {
		if (V <= 0.04045f)
			return V / 12.92f;

		return powf(((V + 0.055f) / 1.055f), 2.4f);
	}

	static float fxyz(float xr) {
		if (xr > epsilonXYZ)
			return powf(xr, 1.0f / 3.0f);

		return 7.787f * xr + 16 / 116.0f;
	}


	void toLAB(float &l, float &a, float &b) const
	{
		float r0 = R / 255.0f;
		float g0 = G / 255.0f;
		float b0 = B / 255.0f;

		r0 = linear(r0);
		g0 = linear(g0);
		b0 = linear(b0);

		float X = 0.4124564f * r0 + 0.3575761f * g0 + 0.1804375f * b0;
		float Y = 0.2126729f * r0 + 0.7151522f * g0 + 0.0721750f * b0;
		float Z = 0.0193339f * r0 + 0.1191920f * g0 + 0.9503041f * b0;

		//Reference white D50
		float Xr = 95.05f / 100.0f;
		float Yr = 100.000f / 100.0f;
		float Zr = 108.88f / 100.0f;

		float xr = X / Xr;
		float yr = Y / Yr;
		float zr = Z / Zr;

		//Lightness
		float L;
		if (yr > epsilonXYZ)
			L = 116 * powf(yr, 1.0f / 3.0f) - 16.0f;
		else
			L = (kXYZ * yr);


		float Al = 500 * (fxyz(xr) - fxyz(yr));
		float B = 200 * (fxyz(yr) - fxyz(zr));

		l = L;
		a = Al;
		b = B;
	}

	double difference(const Color& c2, float s)
	{
		float l1 = 0.0f;
		float a1 = 0.0f;
		float b1 = 0.0f;

		c2.toLAB(l1, a1, b1);

		float l2 = 0.0f;
		float a2 = 0.0f;
		float b2 = 0.0f;

		toLAB(l2, a2, b2);

		double sim = -(powf(l1 - l2, 2.0) + powf(a1 - a2, 2.0) + powf(b1 - b2, 2.0)) / (3.0f * 2.0f * s * s) * 100;
		double dist = expf(sim);
		return dist;
	}

	double e_difference(const Color& c2)
	{
		float l1 = 0.0f;
		float a1 = 0.0f;
		float b1 = 0.0f;

		c2.toLAB(l1, a1, b1);

		float l2 = 0.0f;
		float a2 = 0.0f;
		float b2 = 0.0f;

		toLAB(l2, a2, b2);

		return sqrtf(powf(l1 - l2, 2.0) + powf(a1 - a2, 2.0) + powf(b1 - b2, 2.0));
	}

	float intensity() {
		return (R + G + B) / 3.0f;
	}
};

#endif // ___RGBA_COLOR_IMAGING_H__