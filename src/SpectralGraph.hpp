#ifndef R_1_PUBLIC_GRAPH_SPECTRAL_H_
#define R_1_PUBLIC_GRAPH_SPECTRAL_H_

#include <Eigen/Dense>
#include "Graph.hpp"

#include <algorithm>

namespace R1 {
	class SpectralGraph : public Graph {
	public:
		virtual ~SpectralGraph() {}

		SpectralGraph(int NVertex) : Graph(NVertex) {

		}

		SpectralGraph() {

		}

		//
		// spectral graph properties
		//

		Eigen::VectorXf spectrum() {
			Eigen::SelfAdjointEigenSolver<Eigen::MatrixXf> eg1(adjacency_matrix);
			return eg1.eigenvalues();
		}

		// The energy of G was first defined by Gutman in 1978 as 
		// the sum of the absolute values of its eigenvalues
		double energy() {
			Eigen::MatrixXf::EigenvaluesReturnType eg = adjacency_matrix.eigenvalues();
			double energySum = 0;

			for (int i = 0; i < eg.rows(); i++)	{
				for (int j = 0; j < eg.cols(); j++)	{
					energySum += sqrt(pow(eg(i, j).real(), 2)  + pow(eg(i, j).imag(), 2));
				}
			}

			return energySum;
		}

		// Two graphs are called isospectral or cospectral if the
		// adjacency matrices of the graphs have equal multisets
		// of eigenvalues.
		static bool isospectral(SpectralGraph g1, SpectralGraph g2)
		{
			return g1.spectrum() == g2.spectrum();
		}

		// Two graphs are called isospectral or cospectral if the
		// adjacency matrices of the graphs have equal multisets
		// of eigenvalues.
		static bool isospectral(SpectralGraph g1, SpectralGraph g2, double epsilon)	{
			Eigen::SelfAdjointEigenSolver<Eigen::MatrixXf> eg1(g1.adjacency_matrix);
			Eigen::SelfAdjointEigenSolver<Eigen::MatrixXf> eg2(g2.adjacency_matrix);

			if (eg1.info() != Eigen::Success) {
				return false;
			}

			if (eg2.info() != Eigen::Success) {
				return false;
			}

			Eigen::VectorXf normalizedEg1 = eg1.eigenvalues().normalized();
			Eigen::VectorXf normalizedEg2 = eg2.eigenvalues().normalized();	

			if (normalizedEg1.cols() != normalizedEg2.cols()) {
				return false;
			}

			if (normalizedEg1.rows() != normalizedEg2.rows()) {
				return false;
			}

			for (int i = 0; i < normalizedEg1.rows(); i++)	{
				for (int j = 0; j < normalizedEg1.cols(); j++)	{
					if (abs(normalizedEg1(i, j) - normalizedEg2(i, j)) > epsilon) {
						return false;
					}
				}
			}

			return true;
		}

		//////////////////////////////////////////////////////////////////////////
		//								Laplacian
		//////////////////////////////////////////////////////////////////////////

		float laplacian(int u, int v) {

			float du = degree(u);

			if (u == v && du != 0) {
				return du;
			} else if (adjacency_matrix(u, v) > 0.0f) {
				return -adjacency_matrix(u, v);
			}

			return 0.0f;
		}


		float normalizedLaplacian(int u, int v) {

			float du = degree(u);

			if (u == v && du != 0) {
				return 1.0;
			}
			else if (adjacency_matrix(u, v) > 0.0) {

				float dv = degree(v);

				return -adjacency_matrix(u, v) / sqrt(du * dv);
			}

			return 0.0f;
		}

		Eigen::MatrixXf getLaplacianMatrix() {
			Eigen::MatrixXf degree_matrix(n_vertex, n_vertex);
			degree_matrix.setZero();

			for (int i = 0; i < n_vertex; i++)	{
				degree_matrix(i, i) = degree(i);
			}

			return degree_matrix - adjacency_matrix;
		}

		Eigen::MatrixXf getAjacencyMatrix() {
			return adjacency_matrix;
		}
		

		Eigen::MatrixXf getLaplacianMatrix_obsolete() {
			Eigen::MatrixXf laplacian_matrix(n_vertex, n_vertex);

			for (int i = 0; i < n_vertex; i++)	{
				for (int j = 0; j < n_vertex; j++)	{
					laplacian_matrix(i, j) = laplacian(i, j);
				}

				if (n_vertex > 10 && (i % (n_vertex / 10) == 0))
					std::cout << "LAPLACIAN : " << (100 * i) / n_vertex << std::endl;
			}

			return laplacian_matrix;
		}

		Eigen::MatrixXf getNormalizedLaplacianMatrix() {
			Eigen::MatrixXf degree_matrix(n_vertex, n_vertex);
			degree_matrix.setZero();

			for (int i = 0; i < n_vertex; i++)	{
				degree_matrix(i, i) = 1.0f / (sqrt(degree(i)));
			}
			
			return Eigen::MatrixXf::Identity(n_vertex, n_vertex) - degree_matrix * adjacency_matrix * degree_matrix;
		}

		Eigen::MatrixXf getNormalizedLaplacianMatrix_obsolete() {
			Eigen::MatrixXf laplacian_matrix(n_vertex, n_vertex);

			for (int i = 0; i < n_vertex; i++)	{
				for (int j = 0; j < n_vertex; j++)	{
					laplacian_matrix(i, j) = normalizedLaplacian(i, j);
				}
			}

			return laplacian_matrix;
		}

		Eigen::VectorXcf laplacianSpectrum() {

			Eigen::MatrixXf L = getLaplacianMatrix();
			Eigen::EigenSolver<Eigen::MatrixXf> eg1(L);
			return eg1.eigenvalues();
		}

		Eigen::VectorXcf normalizedLaplacianSpectrum() {

			Eigen::MatrixXf L = getLaplacianMatrix();
			Eigen::EigenSolver<Eigen::MatrixXf> eg1(L);
			return eg1.eigenvalues();
		}

		Eigen::MatrixXcf laplacianEigenVectors() {

			Eigen::MatrixXf L = getLaplacianMatrix();

			Eigen::EigenSolver<Eigen::MatrixXf> eg1(L);
			return eg1.eigenvectors();
		}

		Eigen::MatrixXcf normalizedLaplacianEigenVectors() {

			Eigen::MatrixXf L = getNormalizedLaplacianMatrix();

			Eigen::EigenSolver<Eigen::MatrixXf> eg1(L);
			return eg1.eigenvectors();
		}

		// The Fiedler vector of the graph Laplacian
		Eigen::VectorXf getFiedlerVector(bool normalize = false) {

			Eigen::MatrixXf L;
			if (!normalize)
				L = getLaplacianMatrix();
			else
				L = getNormalizedLaplacianMatrix();

			Eigen::SelfAdjointEigenSolver<Eigen::MatrixXf> eg1;
			eg1.compute(L, Eigen::DecompositionOptions::ComputeEigenvectors);

			Eigen::VectorXf spec = eg1.eigenvalues();

			std::cout << spec << std::endl;

			float m0 = FLT_MAX;
			float m1 = FLT_MAX;

			int m0_index = 0;
			int m1_index = 0;

			for (int i = 0; i < spec.size(); i++)
			{
				if (spec(i) < m0)
				{
					m0 = spec(i);
					m0_index = i;
				}
				else if (spec(i) <= m1) 
				{
					m1 = spec(i);
					m1_index = i;
				}
			}

			return eg1.eigenvectors().col(m1_index);
		}

		void getSortedFiedlerVector(std::vector<float>& fiedlerVector, bool normalize = false)
		{
			Eigen::VectorXf fv = getFiedlerVector(normalize);

			for (int i = 0; i < fv.size(); i++)
			{
				fiedlerVector.push_back(fv(i));
			}
			
			std::sort(fiedlerVector.begin(), fiedlerVector.end());
		}

		static float FiedlerVectorSimilarity(SpectralGraph& sg1, SpectralGraph& sg2)
		{
			std::vector<float> fiedlerVector1;
			std::vector<float> fiedlerVector2;

			sg1.getSortedFiedlerVector(fiedlerVector1);
			sg2.getSortedFiedlerVector(fiedlerVector2);

			size_t n = std::max(fiedlerVector1.size(), fiedlerVector2.size());

			Eigen::VectorXf fv1(n);
			Eigen::VectorXf fv2(n);

			fv1.setZero();
			fv2.setZero();

			for (size_t i = 0; i < n; i++)
			{
				if (i < fiedlerVector1.size())
				{
					fv1(i) = fiedlerVector1[i];
				}

				if (i < fiedlerVector2.size())
				{
					fv2(i) = fiedlerVector2[i];
				}
			}

			return fv1.dot(fv2) / (fv1.norm() * fv2.norm());
		}

		void saveDotFile(std::string filename)
		{
			std::ofstream out(filename);
			out << "graph spectralgraph {\n";
			out << "overlap = scale;\n";

			auto cluster = spectralClusters(false);

			std::unordered_map<int, int> nodeCulster;

			for (size_t i = 0; i < cluster.clusters.size(); i++)
			{
				nodeCulster[cluster.clusters[i].id] = cluster.clusters[i].cluster_id;

				if (cluster.clusters[i].cluster_id == 0)
					out << cluster.clusters[i].id << "[color = blue];" << std::endl;
				else if (cluster.clusters[i].cluster_id == 1)
					out << cluster.clusters[i].id << "[color = green];" << std::endl;
				else 
					out << cluster.clusters[i].id << "[color = red];" << std::endl;
			}

			for (int i = 0; i < n_vertex; i++)
			{
				for (int j = 0; j < i; j++)
				{
					if (is_connected(i, j))
					{
						out << i << "--" << j << ";" << std::endl;
					}
				}
			}

			out << "}\n";

			out.close();
		}

		ClusterCollection spectralClusters(bool normalize = false)
		{
			ClusterCollection cc;

			cc.size = 2;

			Eigen::VectorXf fv = getFiedlerVector(normalize);

			std::cout << fv << std::endl;

			for (int i = 0; i < fv.size(); i++)
			{
				ClusteredNode cnode;
				cnode.id = i;

				if (fabs(fv(i)) < 0.001) {
					cnode.cluster_id = 0;
				} else if (fv(i) > 0.0) {		
					cnode.cluster_id = 1;
				} else {
					cnode.cluster_id = 2;
				}

				cc.clusters.push_back(cnode);
			}

			return cc;
		}

		/*static bool isospectralLaplacian(SpectralGraph g1, SpectralGraph g2, double epsilon) {

			Eigen::VectorXf normalizedEg1 = g1.laplacianSpectrum().normalized();
			Eigen::VectorXf normalizedEg2 = g2.laplacianSpectrum().normalized();	

			if (normalizedEg1.cols() != normalizedEg2.cols()) {
				return false;
			}

			if (normalizedEg1.rows() != normalizedEg2.rows()) {
				return false;
			}

			for (int i = 0; i < normalizedEg1.rows(); i++)	{
				for (int j = 0; j < normalizedEg1.cols(); j++)	{
					if (abs(normalizedEg1(i, j) - normalizedEg2(i, j)) > epsilon) {
						return false;
					}
				}
			}

			return true;
		}

		static double distance(SpectralGraph g1, SpectralGraph g2) {
			Eigen::VectorXf normalizedEg1 = g1.laplacianSpectrum().normalized();
			Eigen::VectorXf normalizedEg2 = g2.laplacianSpectrum().normalized();	

			double dist = 0.0;

			if (normalizedEg1.cols() != normalizedEg2.cols()) {
				return 0.0;
			}

			if (normalizedEg1.rows() != normalizedEg2.rows()) {
				return 0.0;
			}

			for (int i = 0; i < normalizedEg1.rows(); i++)	{
				for (int j = 0; j < normalizedEg1.cols(); j++)	{
					dist += abs(normalizedEg1(i, j) - normalizedEg2(i, j));
				}
			}

			return dist;
		}*/
	};
}

#endif