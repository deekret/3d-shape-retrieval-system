#pragma once

#include <vector>
#include "VectorAttributes.h"
#include "ScalarAttributes.h"
#include <stdio.h>
#include <string>
#include <Eigen/Dense>

using namespace std;


class Grid
{
public:
	Grid(int P, int C): scalars(P), pointNormals(P), faceNormals(C), pointsZ(P)															//Resize the array for the points and cells for the grid
	{
		pointsX.resize(P);
		pointsY.resize(P);
		pointsZ.resize(P);
		cells.resize(3 * C);
	}

	int numPoints()
	{
		return pointsX.size();
	}

	int numCells()
	{
		return cells.size() / 3;
	}

	void getPoint(int i, float* p);

	void setPoint(int i, vector<float> p);

	void setCell(int cell, vector<int> vertices);

	void setClass(string newcls);

	string getClass();

	int	 getCell(int cell, int* vertices);

	int findCell(float* p);

	void normalize();

	void computeFaceNormals();

	void computeVertexNormals();

	void computeCovarianceMatrix();

	void computeEigenvectors();

	VectorAttributes& getFaceNormals()
	{
		return faceNormals;
	}

	VectorAttributes& getPointNormals()
	{
		return pointNormals;
	}

	Eigen::Matrix3f& getCovarianceMatrix() {
		return covarianceMatrix;
	}

	vector<vector<float>> getEigenvectors() {
		vector<vector<float>> vectors;
		vectors.push_back(eigenVec1);
		vectors.push_back(eigenVec2);
		vectors.push_back(eigenVec3);

		return vectors;
	}

	vector<Point3d> getCellCentroids();

	void momentTest();

	int sgn(float x) {
		if (x > 0) return 1;
		if (x < 0) return -1;
		return 1;
	}
	void PCARotation();

	vector<float>		eigenVec1, eigenVec2, eigenVec3;
	int					f0, f1, f2;

protected:

	ScalarAttributes	scalars;

	Eigen::Matrix3f		covarianceMatrix;

	vector<Point3d>		cellCentroids;

	vector<float>		pointsX, pointsY, pointsZ;
	vector<int>			cells;
	VectorAttributes    pointNormals;
	VectorAttributes    faceNormals;
	std::string	cls;

	float				surfaceArea;
	float				compactness;
	float				boundingBoxVolume;
	float				diameter;
	float				eccentricity;

};



