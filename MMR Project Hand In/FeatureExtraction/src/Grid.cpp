#include "..\header\Grid.h"
#include <iostream>
#include <algorithm>
#include <random>

using namespace std;



void Grid::getPoint(int i, float* p)
{
	p[0] = pointsX[i];
	p[1] = pointsY[i];
	p[2] = pointsZ[i];
}


void Grid::setPoint(int i, vector<float> p)
{
	pointsX[i] = p[0];
	pointsY[i] = p[1];
	pointsZ[i] = p[2];
}

int	Grid::getCell(int cell, int* vertices)
{
	vertices[0] = cells[3 * cell];
	vertices[1] = cells[3 * cell + 1];
	vertices[2] = cells[3 * cell + 2];

	/*int len = sizeof(vertices);

	for (int i = 0; i < len; i++) {
		cout << vertices[i] << endl;
		vertices[i] = cells[len * cell + i];
	}*/

	return 3;
}

void Grid::setCell(int cell, vector<int> vertices)
{
	/*cells[3 * cell] = vertices[0];
	cells[3 * cell + 1] = vertices[1];
	cells[3 * cell + 2] = vertices[2];*/

	int len = vertices.size();
	for (int i = 0; i < len; i++) {
		cells[3 * cell + i] = vertices[i];
	}

}

void Grid::normalize()												//Normalize the grid in the [-1,1] cube
{
	float minX = 1.0e6, minY = 1.0e6, minZ = 1.0e6;
	float maxX = -1.0e6, maxY = -1.0e6, maxZ = -1.0e6;

	for (int i = 0; i < numPoints(); ++i)							//1. Determine the bounding-box of all points in the grid
	{
		float p[3];
		getPoint(i, p);
		minX = min(p[0], minX); maxX = max(p[0], maxX);
		minY = min(p[1], minY); maxY = max(p[1], maxY);
		minZ = min(p[2], minZ); maxZ = max(p[2], maxZ);
	}

	float sizeX = maxX - minX;										//2. Compute a single scaling factor that best fits the grid
	sizeX = (sizeX) ? 1 / sizeX : 1;								//   in the [-1,1] cube. Using a single factor for x,y, and z
	float sizeY = maxY - minY;										//   ensures that the object is scaled while keeping its
	sizeY = (sizeY) ? 1 / sizeY : 1;								//   aspect ratio.
	float sizeZ = maxZ - minZ;
	sizeZ = (sizeZ) ? 1 / sizeZ : 1;

	float scale = min(sizeX, min(sizeY, sizeZ));

	for (int i = 0; i < numPoints(); ++i)							//3. Use the scaling factor computed above to scale all grid
	{																//   points in the [-1,1] cube
		float p[3];
		getPoint(i, p);

		p[0] = 2 * ((p[0] - minX) * scale - 0.5);
		p[1] = 2 * ((p[1] - minY) * scale - 0.5);
		p[2] = 2 * ((p[2] - minZ) * scale - 0.5);

		vector<float> v;
		v.push_back(p[0]);
		v.push_back(p[1]);
		v.push_back(p[2]);
		setPoint(i, v);
	}
}

Point3d cellNormal(int size, Point3d* p)						//Compute the normal of a cell whose three vertices are given in p[]
{
	Point3d edge1 = p[1] - p[0];								//We assume that the cell is a triangle. Then, the normal is the
	Point3d edge2 = p[2] - p[1];								//normalized (unit length) value of the cross-product of two cell edges.
	Point3d normal = edge1.cross(edge2);
	normal.normalize();
	return normal;
}

void Grid::computeFaceNormals()								//Compute face normals for the grid. For this, we use the cellNormal()
{															//function presented above, for all grid triangles.
	for (int i = 0; i < numCells(); ++i)
	{
		int cell[10];
		int size = getCell(i, cell);

		Point3d points[10];
		for (int j = 0; j < size; ++j)
		{
			float p[3];
			getPoint(cell[j], p);
			points[j] = Point3d(p);
		}

		Point3d face_normal = cellNormal(size, points);

		faceNormals.setC0Vector(i, face_normal.data);
	}
}

void Grid::computeVertexNormals()							//Compute vertex normals for the grid. For this, we add, to each vertex,
{															//the normals of all cells that use that vertex. Next, we normalize the result.
	for (int i = 0; i < numCells(); ++i)
	{
		int cell[10];
		int size = getCell(i, cell);

		Point3d face_normal = faceNormals.getC0Vector(i);

		for (int j = 0; j < size; ++j)
		{
			Point3d point_normal = pointNormals.getC0Vector(cell[j]);

			point_normal += face_normal;

			pointNormals.setC0Vector(cell[j], point_normal.data);
		}
	}

	for (int i = 0; i < numPoints(); ++i)
	{
		Point3d point_normal = pointNormals.getC0Vector(i);
		point_normal.normalize();
		pointNormals.setC0Vector(i, point_normal.data);
	}
}

void Grid::computeCovarianceMatrix() {

	double means[3] = { 0, 0, 0 };
	vector<vector<float>> covariance;
	vector<vector<float>> points;

	for (int i = 0; i < 3; i++) {
		vector<float> entry;
		for (int j = 0; j < 3; j++) {
			entry.push_back(0.0f);
		}
		covariance.push_back(entry);
	}

	for (int i = 0; i < pointsX.size(); i++) {
		means[0] += pointsX[i];
		means[1] += pointsY[i];
		means[2] += pointsZ[i];
		vector<float> point;

		point.push_back(pointsX[i]);
		point.push_back(pointsY[i]);
		point.push_back(pointsZ[i]);
		points.push_back(point);
	}
		
	means[0] /= points.size(), means[1] /= points.size(), means[2] /= points.size();

	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) {
			covariance[i][j] = 0.0f;
			for (int k = 0; k < points.size(); k++)
				covariance[i][j] += (means[i] - points[k][i]) *
				(means[j] - points[k][j]);
			covariance[i][j] /= points.size() - 1;
			//cout << covariance[i][j];
			//cout << ";";
		}
		//cout << endl;
	}
	covarianceMatrix << covariance[0][0], covariance[0][1], covariance[0][2],
		covariance[1][0], covariance[1][1], covariance[1][2],
		covariance[2][0], covariance[2][1], covariance[2][2];
	//cout << covarianceMatrix << endl;
}

void Grid::computeEigenvectors() {

	computeCovarianceMatrix();

	Eigen::SelfAdjointEigenSolver<Eigen::Matrix3f> eig(covarianceMatrix);
	cout << "Eigenvalues: ";
	cout << eig.eigenvalues() << endl;
	cout << "Eigenvectors: ";
	cout << eig.eigenvectors() << endl << endl;

}

vector<Point3d> Grid::getCellCentroids() {

	vector<Point3d> centroids;
	cellCentroids.clear();

	for (int i = 0; i < numCells(); ++i)
	{
		int cell[10];
		int size = getCell(i, cell);

		Point3d points[10];
		Point3d centroid;
		for (int j = 0; j < size; ++j)
		{
			float p[3];
			getPoint(cell[j], p);
			centroid.x += Point3d(p).x;
			centroid.y += Point3d(p).y;
			centroid.z += Point3d(p).z;
		}
		centroid.x /= size;
		centroid.y /= size;
		centroid.z /= size;
		cellCentroids.push_back(centroid);
		centroids.push_back(centroid);
	}

	return centroids;
}

void Grid::momentTest() {

	Eigen::Matrix3i F;

	getCellCentroids();

	float fX = 0;
	float fY = 0;
	float fZ = 0;

	for (int i = 0; i < numCells(); i++) {

		fX += (sgn(cellCentroids[i].x) * pow(cellCentroids[i].x, 2.0));
		fY += (sgn(cellCentroids[i].y) * pow(cellCentroids[i].y, 2.0));
		fZ += (sgn(cellCentroids[i].z) * pow(cellCentroids[i].z, 2.0));
	}

	cout << sgn(fX) << endl;
	cout << sgn(fY) << endl;
	cout << sgn(fZ) << endl;

	for (int i = 0; i < numPoints(); i++) {
		pointsX[i] = sgn(fX) * pointsX[i];
		pointsY[i] = sgn(fY) * pointsY[i];
		pointsZ[i] = sgn(fZ) * pointsZ[i];
	}

	F << sgn(fX), 0, 0,
		0, sgn(fY), 0,
		0, 0, sgn(fZ);

	cout << sgn(fX) << endl;
	cout << sgn(fY) << endl;
	cout << sgn(fZ) << endl;

}

void Grid::PCARotation() {
	computeCovarianceMatrix();
	computeEigenvectors();
	Eigen::SelfAdjointEigenSolver<Eigen::Matrix3f> eig(covarianceMatrix);
	vector<float> vec1, vec2, vec3;

	for (int i = 0; i < 3; i++) {
		vec1.push_back(eig.eigenvectors().col(0)(i));
		vec2.push_back(eig.eigenvectors().col(1)(i));
		vec3.push_back(eig.eigenvectors().col(2)(i));
	}


	for (int i = 0; i < numPoints(); i++) {

		Point3d newPoint;

		newPoint.x = (vec1[0] * pointsX[i] + vec1[1] * pointsY[i] + vec1[2] * pointsZ[i]);
		newPoint.y = (vec2[0] * pointsX[i] + vec2[1] * pointsY[i] + vec2[2] * pointsZ[i]);


		//vec3 = crossproduct of vec1 and vec2
		/*vec1[0] = vec3[1] * vec2[2] - vec3[2] * vec2[1];
		vec1[0] = vec3[2] * vec2[0] - vec3[0] * vec2[2];
		vec1[0] = vec3[0] * vec2[1] - vec3[1] * vec2[0];*/

		newPoint.x = (vec3[0]	* pointsX[i]	+ vec3[1]	* pointsY[i]	+ vec3[2]	* pointsZ[i]);
		newPoint.y = (vec2[0]	* pointsX[i]	+ vec2[1]	* pointsY[i]	+ vec2[2]	* pointsZ[i]);
		newPoint.z = (vec1[0]	* pointsX[i]	+ vec1[1]	* pointsY[i]	+ vec1[2]	* pointsZ[i]);

		pointsX[i] = newPoint.x;
		pointsY[i] = newPoint.y;
		pointsZ[i] = newPoint.z;
	}


}


string Grid::getClass()
{
	return cls;
}

void Grid::setClass(string newcls)
{
	cls = newcls;
}

float Grid::calculateSurfaceArea() {

	float area = 0;

	for (int i = 0; i < numCells(); i++) {
		
		float cellArea = 0;
		float s = 0;
		Point3d points[3];

		int vertices[3];
		int size = getCell(i, vertices);

		if (size != 3) {
			cout << "NOT A TRIANGLE CELL! ABORT!";
			break;
		}

		float d1, d2, d3;
		float p0[3];
		float p1[3];
		float p2[3];

		getPoint(vertices[0], p0);
		getPoint(vertices[1], p1);
		getPoint(vertices[2], p2);
		/*cout << Point3d(p0).x << " " << Point3d(p0).y << " " << Point3d(p0).z << endl;
		cout << Point3d(p1).x << " " << Point3d(p1).y << " " << Point3d(p1).z << endl;
		cout << Point3d(p2).x << " " << Point3d(p2).y << " " << Point3d(p2).z << endl;*/


		d1 = sqrt(pow((Point3d(p1).x - Point3d(p0).x),2.0) + pow((Point3d(p1).y - Point3d(p0).y),2.0) + pow((Point3d(p1).z - Point3d(p0).z),2.0));
		d2 = sqrt(pow((Point3d(p2).x - Point3d(p1).x),2.0) + pow((Point3d(p2).y - Point3d(p1).y),2.0) + pow((Point3d(p2).z - Point3d(p1).z),2.0));
		d3 = sqrt(pow((Point3d(p0).x - Point3d(p2).x),2.0) + pow((Point3d(p0).y - Point3d(p2).y),2.0) + pow((Point3d(p0).z - Point3d(p2).z),2.0));
		s = (d1 + d2 + d3) / 2;

		cellArea = sqrt(s * ((s - d1) * (s - d2) * (s - d3)));

		area += cellArea;
	}
	surfaceArea = area;
	return area;
}

float Grid::calculateVolume() {

	float tempVolume = 0;

	for (int i = 0; i < numCells(); i++) {

		float pv = 0;
		float s = 0;
		Point3d points[3];

		int vertices[3];
		int size = getCell(i, vertices);

		if (size != 3) {
			cout << "NOT A TRIANGLE CELL! ABORT!";
			break;
		}

		float d1, d2, d3;
		float p0[3];
		float p1[3];
		float p2[3];

		getPoint(vertices[0], p0);
		getPoint(vertices[1], p1);
		getPoint(vertices[2], p2);

		points[0] = p0;
		points[1] = p1;
		points[2] = p2;

		/*cout << cellNormal(i, points).x << endl;
		cout << cellNormal(i, points).y << endl;
		cout << cellNormal(i, points).z << endl;*/

		pv = Point3d(p0).x * Point3d(p1).y * Point3d(p2).z
			+ Point3d(p0).y * Point3d(p1).z * Point3d(p2).x
			+ Point3d(p0).z * Point3d(p1).x * Point3d(p2).y
			- Point3d(p0).x * Point3d(p1).z * Point3d(p2).y
			- Point3d(p0).y * Point3d(p1).x * Point3d(p2).z
			- Point3d(p0).z * Point3d(p1).y * Point3d(p2).x;

		tempVolume += pv;
	}

	volume = tempVolume / 6;
	if (sgn(volume) == -1) {
		volume *= -1;
	}
	return volume;
}

float Grid::calculateSphericity() {

	calculateSurfaceArea();
	calculateVolume();

	compactness = pow(surfaceArea, 3.0) / (36 * M_PI * pow(volume, 2.0));
	sphericity = 1 / compactness;
	return compactness ;
}

float Grid::calculateBoundingBoxVol() {

	float s1, s2, s3;

	s1 = (maxX - minX);
	s2 = (maxY - minY);
	s3 = (maxZ - minZ);

	boundingBoxVolume = s1 * s2 * s3;
	return boundingBoxVolume;
}

float Grid::calculateDiameter() {

	float disFromBary = -1;
	float disFromFarthest = -1;
	float far_x = 0;
	float far_y = 0;
	float far_z = 0;

	for (int i = 0; i < numPoints(); i++) {				// Find the point furthest from the barycenter

		float p[3];
		getPoint(i, p);

		float d = 0;
		d = sqrt(pow(p[0], 2.0) + pow(p[1], 2.0) + pow(p[2], 2.0)); //Distance From barycenter
		if (d > disFromBary)
		{
			disFromBary = d;
			far_x = p[0];
			far_y = p[1];
			far_z = p[2];
		}
	}

	for (int i = 0; i < numPoints(); i++) {				// Find the point furthest from the previously found point

		float p[3];
		getPoint(i, p);

		float d = 0;

		d = sqrt(pow((p[0] - far_x), 2.0) + pow((p[1] - far_y), 2.0) + pow((p[2] - far_z), 2.0));

		if (d > disFromFarthest)
		{
			disFromFarthest = d;
		}

	}
	diameter = disFromFarthest;
	return disFromFarthest;
}

float Grid::calculateEccentricity() {
	float majorEigenValue;
	float minorEigenValue;
	float ecc;
	computeCovarianceMatrix();

	Eigen::SelfAdjointEigenSolver<Eigen::Matrix3f> eig(covarianceMatrix);


	majorEigenValue = eig.eigenvalues()[0];
	minorEigenValue = eig.eigenvalues()[0];

	for (int i = 1; i < 3; i++) {
		if (majorEigenValue < eig.eigenvalues()[i])
		{
			majorEigenValue = eig.eigenvalues()[i];

		}
		if (minorEigenValue > eig.eigenvalues()[i])
		{
			minorEigenValue = eig.eigenvalues()[i];
		}
	}


	//cout << abs(majorEigenValue)/abs(minorEigenValue) << endl;
	ecc = abs(majorEigenValue) / abs(minorEigenValue);
	eccentricity = ecc;

	return ecc;
}

void Grid::calculateA3(int n) 
{
	double res;
	float angle;
	int k = pow(n, 1.0 / 3.0);

	for (int i = 0; i < k; i++) {

		int p1 = rand() % numPoints();

		for (int j = 0; j < k; j++) {

			int p2 = rand() % numPoints();
			if (p1 == p2) {
				continue;											// do not allow equal points;
			}

			for (int l = 0; l < k; l++) {

				int p3 = rand() % numPoints();

				if (p1 == p3 || p2 == p3)
				{
					continue;										// do not allow equal points;
				}

				double ba[3] = { pointsX[p1] - pointsX[p2], pointsY[p1] - pointsY[p2], pointsZ[p1] - pointsZ[p2] };
				
				double bc[3] = { pointsX[p3] - pointsX[p2], pointsY[p3] - pointsY[p2], pointsZ[p3] - pointsZ[p2] };
				
				//normalizing
				double baVec = sqrt(ba[0] * ba[0] + ba[1] * ba[1] + ba[2] * ba[2]);
				double bcVec = sqrt(bc[0] * bc[0] + bc[1] * bc[1] + bc[2] * bc[2]);
				
				double baNorm[3] = { ba[0] / baVec, ba[1] / baVec, ba[2] / baVec };
				double bcNorm[3] = { bc[0] / bcVec, bc[1] / bcVec, bc[2] / bcVec };
				
				//calculating the dot product
				res = baNorm[0] * bcNorm[0] + baNorm[1] * bcNorm[1] + baNorm[2] * bcNorm[2];
				angle = acos(res) * 180.0 / M_PI;

				if (angle < A3min)
					A3min = angle;
				if (angle > A3max)
					A3max = angle;
				A3s.push_back(angle);
			}
		}
	}
}

void Grid::calculateD1() {

	// Calculates the distance between the barycenter and a random vertex.

	float distance;
	D1max = FLT_MIN;
	D1min = FLT_MAX;

	for (int i = 0; i < numPoints(); i++) 
	{
		distance = sqrt(pow(pointsX[i], 2.0) + pow(pointsY[i], 2.0) + pow(pointsZ[i], 2.0));
		if (distance < D1min)
			D1min = distance;
		if (distance > D1max)
			D1max = distance;
		D1s.push_back(distance);
	}
}

void Grid::calculateD2(int n) {

	// Calculate distance between two random vertices.
	// The distances are added into the histogram vector D2hist with 10 bins, each bin being 0.15.
	// Max distance = 1.75?

	float distance;
	D1max = FLT_MIN;
	D1min = FLT_MAX;

	int k = pow(n, 1.0 / 2.0);

	for (int i = 0; i < k; i++) {

		int p1 = rand() % numPoints();

		for (int j = 0; j < k; j++) {

			int p2 = rand() % numPoints();

			if (p1 == p2) {
				continue;										// do not allow equal points;
			}

			distance = sqrt(pow((pointsX[p1] - pointsX[p2]), 2.0) + pow((pointsY[p1] - pointsY[p2]), 2.0)
				+ pow((pointsZ[p1] - pointsZ[p2]), 2.0));
			D2s.push_back(distance);

			if (distance < D2min)
				D2min = distance;
			if (distance > D2max)
				D2max = distance;
		}
	}
}

void Grid::calculateD3(int n) {

	// Calculate distance between two random vertices.
	// The areas are added into the histogram vector D2hist with 10 bins, each bin being 0.15.
	// Max area = 0.85?

	float area, s, d1, d2, d3;

	int k = pow(n, 1.0 / 3.0);

	for (int i = 0; i < k; i++) {

		int p1 = rand() % numPoints();

		for (int j = 0; j < k; j++) {

			int p2 = rand() % numPoints();
			if (p1 == p2) {
				continue;											// do not allow equal points;
			}

			for (int l = 0; l < k; l++) {

				int p3 = rand() % numPoints();
				if (p1 == p3 || p2 == p3) {
					continue;										// do not allow equal points;
				}

				d1 = sqrt(pow((pointsX[p2] - pointsX[p1]),2.0) + pow((pointsY[p2] - pointsY[p1]),2.0) + pow((pointsZ[p2] - pointsZ[p1]),2.0));
				d2 = sqrt(pow((pointsX[p3] - pointsX[p2]),2.0) + pow((pointsY[p3] - pointsY[p2]),2.0) + pow((pointsZ[p3] - pointsZ[p2]),2.0));
				d3 = sqrt(pow((pointsX[p1] - pointsX[p3]),2.0) + pow((pointsY[p1] - pointsY[p3]),2.0) + pow((pointsZ[p1] - pointsZ[p3]),2.0));
				s = (d1 + d2 + d3) / 2;

				area = sqrt(s * ((s - d1) * (s - d2) * (s - d3)));
				area = sqrt(area);

				D3s.push_back(area);

				if (area < D3min)
					D3min = area;
				if (area > D3max)
					D3max = area;
			}
		}
	}
}

void Grid::calculateD4(int n) {

	// Calculate volume between four random vertices.
	// The volumes are added into the histogram vector D4hist with 10 bins, each bin being 0.033.
	// Max volume = 0.33?

	float volume;
	int k = pow(n, 1.0 / 4.0);

	for (int i = 0; i < k; i++) {

		int p1 = rand() % numPoints();

		for (int j = 0; j < k; j++) {

			int p2 = rand() % numPoints();

			if (p1 == p2) {
				continue;										// do not allow equal points;
			}

			for (int l = 0; l < k; l++)
			{
				int p3 = rand() % numPoints();

				if (p3 == p1 || p3 == p2)
					continue;

				for (int m = 0; m < k; m++)
				{
					int p4 = rand() % numPoints();

					if (p4 == p1 || p4 == p2 || p4 == p3)
						continue;

					float* temp = new float[3];
					getPoint(p1, temp);
					Point3d P1 = Point3d((const float*)temp);
					getPoint(p2, temp);
					Point3d P2 = Point3d((const float*)temp);
					getPoint(p3, temp);
					Point3d P3 = Point3d((const float*)temp);
					getPoint(p4, temp);
					Point3d P4 = Point3d((const float*)temp);
					delete[] temp;

					Point3d numerator = (P1 - P4).dot((P2 - P4).cross(P3 - P4));
					volume = numerator.norm() / 6.0f;
					volume = pow(volume, 1.0 / 3.0);

					D4s.push_back(volume);

					if (volume < D4min)
						D4min = volume;
					if (volume > D4max)
						D4max = volume;
				}
			}
		}
	}
}

vector<vector<float>> Grid::calcHists(float D1min, float D1max, float D2min, float D2max, float D3min, float D3max, float D4min, float D4max, float A3min, float A3max, int bins, string csv)
{
	std::ifstream fin(csv);
	vector<vector<float>> hists;
	for (int i = 0; i < 5; i++)
	{
		vector<float> temp;
		for (int j = 0; j < bins; j++)
			temp.push_back(0);
		hists.push_back(temp);
	}

	vector<float> mins, maxs;
	mins.push_back(D1min);
	mins.push_back(D2min);
	mins.push_back(D3min);
	mins.push_back(D4min);
	mins.push_back(A3min);
	maxs.push_back(D1max);
	maxs.push_back(D2max);
	maxs.push_back(D3max);
	maxs.push_back(D4max);
	maxs.push_back(A3max);


	string line;
	getline(fin, line);

	vector<float> binSizes;
	for (int i = 0; i < 5; i++)
		binSizes.push_back((maxs[i] - mins[i]) / bins);

	vector<int> counts;
	for (int i = 0; i < 5; i++)
		counts.push_back(0);

	while (!fin.eof())
	{
		getline(fin, line);

		vector<string> splits;
		string item;

		stringstream ss(line);
		while (getline(ss, item, ';'))
			splits.push_back(item);

		for (int i = 0; i < 5; i++)
		{
			if (splits.size() <= i)
				break;
			if (splits[i] == "")
				continue;

			float val = stof(splits[i]);
			counts[i]++;

			for (int b = 0; b < bins; b++)
			{
				if (val >= (b * binSizes[i] + mins[i]) && val < (b + 1) * binSizes[i] + mins[i])
				{
					hists[i][b]++;
					break;
				}
				else if (val >= bins * binSizes[i] + mins[i])
				{
					hists[i][bins - 1];
					break;
				}
			}
		}
	}

	for (int i = 0; i < 5; i++)
	{
		for (int j = 0; j < bins; j++)
		{
			hists[i][j] /= counts[i];
		}
	}

	return hists;
}

void Grid::outputHist(string name, int count)
{
	ofstream fout;
	fout.open(name);

	fout << "sep=;" << endl;;

	for (int i = 0; i < count; i++)
	{
		if (D1s.size() > i)
			fout << D1s[i];
		fout << ";";

		if (D2s.size() > i)
			fout << D2s[i];
		fout << ";";

		if (D3s.size() > i)
			fout << D3s[i];
		fout << ";";

		if (D4s.size() > i)
			fout << D4s[i];
		fout << ";";

		if (A3s.size() > i)
			fout << A3s[i];
		fout << endl;
	}

	fout.close();
}



