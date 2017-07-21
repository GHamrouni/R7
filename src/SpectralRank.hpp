#pragma once

#include "Eigen/Dense"

using namespace std;

namespace SpectralRanking
{
	using namespace Eigen;

	MatrixXf SetStochastic(MatrixXf S)
	{
		for (int i = 0; i < S.rows(); i++)
		{
			float sum = S.row(i).sum();

			assert(sum == sum);

			if (sum != 0)
				S.row(i) = S.row(i) / sum;
		}

		return S;
	}

	VectorXf dominantEigenVector(MatrixXf G, float epsilon, int maxRec)
	{
		VectorXf v(G.rows());
		v.fill(1);

		MatrixXf  vT = v.transpose();
		MatrixXf  vT_old(vT);

		for (int i = 0; i<maxRec;i++)
		{
			vT = vT * G;

			if (vT_old.isApprox(vT, epsilon))
				break;

			vT_old = vT;
		}

		v = vT.transpose();

		float norm = v.sum();

		if (norm != 0.0f)
			v = v/norm;

		return v;
	}

	VectorXf PageRank(MatrixXf S, float dampingFactor, float epsilon, int maxRec)
	{
		S = SetStochastic(S);
		MatrixXf I(S.rows(), S.cols());

		I.fill(1.0f/((float)S.rows()));

		MatrixXf G =  (1 - dampingFactor) * I + dampingFactor * S;

		return dominantEigenVector(G, epsilon, maxRec);
	}

	VectorXf CheiRank(MatrixXf S, float dampingFactor, float epsilon, int maxRec)
	{
		S.transposeInPlace();

		return PageRank(S, dampingFactor, epsilon, maxRec);
	}
}
