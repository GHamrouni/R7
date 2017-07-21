/*

Copyright (c) 2013  Ghassen Hamrouni

Abstract:

This module provides an implementation of an Histogram

Author:

Ghassen Hamrouni <ghamrouni.iptech@gmail.com> 23-07-2013

Revision History:

*/

#ifndef R_1_PUBLIC_HISTOGRAM_H_
#define R_1_PUBLIC_HISTOGRAM_H_

#include <algorithm>
#include <vector>

namespace R1 {
	class Histogram
	{
		size_t BinSize;

		float Min;

		float Max;

		std::vector<float> buffer;

		Histogram(size_t BinSize, float min, float max)
		{
			this->BinSize = BinSize;

			this->Max = max;
			this->Min = min;

			//buffer = new float[BinSize];

			for (size_t i = 0; i < BinSize; i++)
			{
				buffer.push_back(0.0);
			}
		}

	public:

		~Histogram()
		{
			//delete buffer;
		}

		static Histogram BuildHistogram(int BinSize, float min, float max)
		{
			Histogram histo = Histogram(BinSize, min, max);
			return histo;
		}

		static Histogram* BuildHistogram(int BinSize, std::vector<float>& data)
		{
			float min = 0;
			float max = 0;

			if (data.size() > 0)
			{
				min = *(std::min_element(data.begin(), data.end()));
				max = *(std::max_element(data.begin(), data.end()));
			}

			Histogram* histo = new Histogram(BinSize, min, max);

			for (size_t i = 0; i < data.size(); i++)
			{
				histo->Add(data[i]);
			}

			return histo;
		}

		void Log10()
		{
			for (size_t i = 0; i < BinSize; i++)
			{
				if (buffer[i] > 0.0)
					buffer[i] = log10(buffer[i]);
			}
		}

		void Log()
		{
			for (size_t i = 0; i < BinSize; i++)
			{
				if (buffer[i] > 0.0)
					buffer[i] = log(buffer[i]);
			}
		}

		float Average() const
		{
			float sum = 0.0;

			for (size_t i = 0; i < BinSize; i++)
			{
				sum += buffer[i];
			}

			sum = sum / (float)BinSize;

			return sum;
		}

		// Normalize with the sum
		void NormalizeSum()
		{
			float sum = 0.0;

			for (size_t i = 0; i < BinSize; i++)
			{
				sum += buffer[i];
			}

			if (sum != 0.0)
			{
				for (size_t i = 0; i < BinSize; i++)
				{
					buffer[i] = buffer[i] / sum;
				}
			}
		}

		/// <summary>
		/// Normalize with the mean
		/// </summary>
		void Normalize()
		{
			float sum = 0.0;

			for (size_t i = 0; i < BinSize; i++)
			{
				sum += buffer[i];
			}

			sum = sum / (float)BinSize;

			if (sum != 0.0)
			{
				for (size_t i = 0; i < BinSize; i++)
				{
					buffer[i] = buffer[i] / sum;
				}
			}
		}

		// Normalize with the max
		void NormalizeMax()
		{
			float sum = *(std::max_element(buffer.begin(), buffer.end()));

			if (sum != 0.0)
			{
				for (size_t i = 0; i < BinSize; i++)
				{
					buffer[i] = buffer[i] / sum;
				}
			}
		}

		// Chi-square correlation coefficient
		float DistanceChiSquare(Histogram* hist)
		{
			float xm = Average();
			float ym = hist->Average();

			float nx = 0.0;

			for (size_t i = 0; i < buffer.size(); i++)
			{
				if ((buffer[i] + hist->buffer[i]) != 0.0)
				{
					nx += pow(buffer[i] - hist->buffer[i], 2.0) / (buffer[i] + hist->buffer[i]);
				}
			}

			return nx * 0.5;
		}

		// Pearson product-moment correlation coefficient
		float DistanceCor(Histogram* hist)
		{
			float xm = Average();
			float ym = hist->Average();

			float s = 0.0f;
			float nx = 0.0f;
			float ny = 0.0f;

			for (size_t i = 0; i < buffer.size(); i++)
			{
				s += (buffer[i] - xm) * (hist->buffer[i] - ym);
				nx += powf(buffer[i] - xm, 2.0f);
				ny += powf(hist->buffer[i] - ym, 2.0f);
			}

			if (nx * ny == 0.0f && (nx != 0.0f || ny != 0.0f))
			{
				return 1.0f;
			}

			return 1.0f - s / sqrtf(nx * ny);
		}

		void Add(float value)
		{
			float v = value - Min;
			float I = Max - Min;

			if (I == 0.0)
			{
				return;
			}

			float index = (v / I) * (BinSize - 1);

			buffer[(int)index] += 1.0;
		}

		void Add(float index_value, float value)
		{
			float v = index_value - Min;
			float I = Max - Min;

			if (I == 0.0)
			{
				return;
			}

			float index = (v / I) * (BinSize - 1);

			buffer[(int)index] += value;
		}
	};
}

#endif // R_1_PUBLIC_HISTOGRAM_H_