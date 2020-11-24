#pragma once

#include <vector>
#include "VectorAttributes.h"
#include "ScalarAttributes.h"
#include <stdio.h>
#include <string>
#include <Eigen/Dense>
#include <fstream>
#include <iostream>

using namespace std;
#ifndef M_PI
#define M_PI (3.14159265358979323846)
#endif

class Grid
{
public:
	Grid(int P, int C) : scalars(P), pointNormals(P), faceNormals(C), pointsZ(P)															//Resize the array for the points and cells for the grid
	{
		pointsX.resize(P);
		pointsY.resize(P);
		pointsZ.resize(P);
		cells.resize(3 * C);

		for (int i = 0; i < 14; i++) {
			D1hist.push_back(0);
			D2hist.push_back(0);
			D3hist.push_back(0);
			D4hist.push_back(0);
			A3hist.push_back(0);
		}
	}

	float calculateSurfaceArea();
	float calculateEccentricity();

	float calculateBoundingBoxVol();
	float calculateDiameter();

	float calculateVolume();

	float calculateSphericity();
	double calculateAngleBetweenPoints();

	float D1min, D1max, D2min, D2max, D3min, D3max, D4min, D4max, A3min, A3max;
	void outputHist(string name, int count);
	static vector<vector<float>> calcHists(float D1min, float D1max, float D2min, float D2max, float D3min, float D3max, float D4min, float D4max, float A3min, float A3max, int bins, string csv);


	void setExtremes(float mix, float max, float miy, float may, float miz, float maz) {
		minX = mix;
		maxX = max;
		minY = miy;
		maxY = may;
		minZ = miz;
		maxZ = maz;
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

	void calculateD1();
	vector<float> getD1hist(float min, float max, int bins) 
	{
		float binSize = (max - min) / bins;

		for (int i = 0; i < D1s.size(); i++)
		{
			float distance = D1s[i];
			for (int b = 0; b < bins; b++) {
				if (distance >= (b * binSize + min) && distance < (b + 1) * binSize + min) {
					D1hist[b] += 1;
					break;
				}
				else if (distance >= bins * binSize + min) {
					D1hist[bins - 1] += 1; 
					break;
				}
			}
		}

		for (int i = 0; i < bins; i++) {
			D1hist[i] = D1hist[i] / D1s.size();
		}
		
		return D1hist;
	}

	vector<float> getD2hist(float min, float max, int bins)
	{
		float binSize = (max - min) / bins;

		for (int i = 0; i < D2s.size(); i++)
		{
			float distance = D2s[i];
			for (int b = 0; b < bins; b++) {
				if (distance >= (b * binSize + min) && distance < (b + 1) * binSize + min) {
					D2hist[b] += 1;
					break;
				}
				else if (distance >= bins * binSize + min) {
					D2hist[bins - 1] += 1;
					break;
				}
			}
		}

		for (int i = 0; i < bins; i++) {
			D2hist[i] = D2hist[i] / D2s.size();
		}

		return D2hist;
	}

	vector<float> getD3hist(float min, float max, int bins)
	{
		float binSize = (max - min) / bins;

		for (int i = 0; i < D3s.size(); i++)
		{
			float distance = D3s[i];
			for (int b = 0; b < bins; b++) {
				if (distance >= (b * binSize + min) && distance < (b + 1) * binSize + min) {
					D3hist[b] += 1;
					break;
				}
				else if (distance >= bins * binSize + min) {
					D3hist[bins - 1] += 1;
					break;
				}
			}
		}

		for (int i = 0; i < bins; i++) {
			D3hist[i] = D3hist[i] / D3s.size();
		}

		return D3hist;
	}

	vector<float> getD4hist(float min, float max, int bins)
	{
		float binSize = (max - min) / bins;

		for (int i = 0; i < D4s.size(); i++)
		{
			float distance = D4s[i];
			for (int b = 0; b < bins; b++) {
				if (distance >= (b * binSize + min) && distance < (b + 1) * binSize + min) {
					D4hist[b] += 1;
					break;
				}
				else if (distance >= bins * binSize + min) {
					D4hist[bins - 1] += 1;
					break;
				}
			}
		}

		for (int i = 0; i < bins; i++) {
			D4hist[i] = D4hist[i] / D4s.size();
		}

		return D4hist;
	}

	vector<float> getA3hist(float min, float max, int bins)
	{
		float binSize = (max - min) / bins;

		for (int i = 0; i < A3s.size(); i++)
		{
			float distance = A3s[i];
			for (int b = 0; b < bins; b++) {
				if (distance >= (b * binSize + min) && distance < (b + 1) * binSize + min) {
					A3hist[b] += 1;
					break;
				}
				else if (distance >= bins * binSize + min) {
					A3hist[bins - 1] += 1;
					break;
				}
			}
		}

		for (int i = 0; i < bins; i++) {
			A3hist[i] = A3hist[i] / A3s.size();
		}

		return A3hist;
	}

	static vector<float> getD1hist(float min, float max, int bins, string csv)
	{
		std::ifstream fin(csv);
		vector<float> D1hist;
		for (int i = 0; i < bins; i++)
			D1hist.push_back(0);

		string line;
		getline(fin, line);

		float binSize = (max - min) / bins;
		int count;

		for (count = 0; !fin.eof(); count++)
		{
			getline(fin, line);

			vector<string> splits;
			string item;

			stringstream ss(line);
			while (getline(ss, item, ';'))
				splits.push_back(item);

			if (splits.size() == 0 || splits[0] == "")
				break;
			float distance = stof(splits[0]);
			for (int b = 0; b < bins; b++) {
				if (distance >= (b * binSize + min) && distance < (b + 1) * binSize + min) {
					D1hist[b] += 1;
					break;
				}
				else if (distance >= bins * binSize + min) {
					D1hist[bins - 1] += 1;
					break;
				}
			}
		}

		for (int i = 0; i < bins; i++) {
			D1hist[i] = D1hist[i] / count;
		}

		return D1hist;
	}

	void calculateD2(int n);
	static vector<float> getD2hist(float min, float max, int bins, string csv) 
	{
		std::ifstream fin(csv);
		vector<float> D2hist;
		for (int i = 0; i < bins; i++)
			D2hist.push_back(0);

		string line;
		getline(fin, line);

		float binSize = (max - min) / bins;
		int count;

		for (count = 0; !fin.eof(); count++)
		{
			getline(fin, line);

			vector<string> splits;
			string item;

			stringstream ss(line);
			while (getline(ss, item, ';'))
				splits.push_back(item);

			if (splits.size() <= 1)
				break;

			string temp = splits[1];
			if (temp == "")
				break;
			float distance = stof(temp);
			for (int b = 0; b < bins; b++) {
				if (distance >= (b * binSize + min) && distance < (b + 1) * binSize + min) {
					D2hist[b] += 1;
					break;
				}
				else if (distance >= bins * binSize + min) {
					D2hist[bins - 1] += 1;
					break;
				}
			}
		}

		for (int i = 0; i < bins; i++) {
			D2hist[i] = D2hist[i] / count;
		}

		return D2hist;
	}

	void calculateD3(int n);
	static vector<float> getD3hist(float min, float max, int bins, string csv) 
	{
		ifstream fin(csv);
		vector<float> D3hist;
		for (int i = 0; i < bins; i++)
			D3hist.push_back(0);

		string line;
		getline(fin, line);

		float binSize = (max - min) / bins;
		int count;

		for (count = 0; !fin.eof(); count++)
		{
			getline(fin, line);

			vector<string> splits;
			string item;

			stringstream ss(line);
			while (getline(ss, item, ';'))
				splits.push_back(item);

			if (splits.size() <= 2)
				break;

			string temp = splits[2];
			if (temp == "")
				break;

			float distance = stof(temp);
			for (int b = 0; b < bins; b++) {
				if (distance >= (b * binSize + min) && distance < (b + 1) * binSize + min) {
					D3hist[b] += 1;
					break;
				}
				else if (distance >= bins * binSize + min) {
					D3hist[bins - 1] += 1;
					break;
				}
			}
		}

		for (int i = 0; i < bins; i++) {
			D3hist[i] = D3hist[i] / count;
		}

		return D3hist;
	}

	void calculateD4(int n);
	static vector<float> getD4hist(float min, float max, int bins, string csv) 
	{
		ifstream fin(csv);
		vector<float> D4hist;
		for (int i = 0; i < bins; i++)
			D4hist.push_back(0);

		string line;
		getline(fin, line);

		float binSize = (max - min) / bins;
		int count;

		for (count = 0; !fin.eof(); count++)
		{
			getline(fin, line);

			vector<string> splits;
			string item;

			stringstream ss(line);
			while (getline(ss, item, ';'))
				splits.push_back(item);

			if (splits.size() <= 3)
				break;

			string temp = splits[3];
			if (temp == "")
				break;

			float distance = stof(temp);
			for (int b = 0; b < bins; b++) {
				if (distance >= (b * binSize + min) && distance < (b + 1) * binSize + min) {
					D4hist[b] += 1;
					break;
				}
				else if (distance >= bins * binSize + min) {
					D4hist[bins - 1] += 1;
					break;
				}
			}
		}

		for (int i = 0; i < bins; i++) {
			D4hist[i] = D4hist[i] / count;
		}

		return D4hist;
	}

	void calculateA3(int n);
	static vector<float> getA3hist(float min, float max, int bins, string csv) 
	{
		ifstream fin(csv);
		vector<float> A3hist;
		for (int i = 0; i < bins; i++)
			A3hist.push_back(0);

		string line;
		getline(fin, line);

		float binSize = (max - min) / bins;
		int count;

		for (count = 0; !fin.eof(); count++)
		{
			getline(fin, line);

			vector<string> splits;
			string item;

			stringstream ss(line);
			while (getline(ss, item, ';'))
				splits.push_back(item);

			if (splits.size() <= 4)
				break;

			string temp = splits[4];
			if (temp == "")
				break;

			float distance = stof(temp);
			for (int b = 0; b < bins; b++) {
				if (distance >= (b * binSize + min) && distance < (b + 1) * binSize + min) {
					A3hist[b] += 1;
					break;
				}
				else if (distance >= bins * binSize + min) {
					A3hist[bins - 1] += 1;
					break;
				}
			}
		}

		for (int i = 0; i < bins; i++) {
			A3hist[i] = A3hist[i] / count;
		}

		return A3hist;
	}

	float				surfaceArea;
	float				volume;
	float				compactness;
	float				sphericity;
	float				minX, maxX, minY, maxY, minZ, maxZ;
	float				boundingBoxVolume;
	float				diameter;
	float				eccentricity;
	vector<float>		D1hist, D2hist, D3hist, D4hist, A3hist;
	string				name, shape;


protected:

	ScalarAttributes	scalars;

	Eigen::Matrix3f		covarianceMatrix;

	vector<float>		eigenVec1, eigenVec2, eigenVec3;

	vector<Point3d>		cellCentroids;

	vector<float>		pointsX, pointsY, pointsZ;
	vector<int>			cells;
	VectorAttributes    pointNormals;
	VectorAttributes    faceNormals;
	std::string	cls;

	vector<float>		D1s, D2s, D3s, D4s, A3s;
};



