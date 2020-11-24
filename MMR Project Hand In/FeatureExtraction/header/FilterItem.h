#pragma once

#include <string>

struct FilterItem
{
	std::string cls;
	int numFaces;
	int numVertices;
	std::string typeOfFace;
	float minX, maxX, minY, maxY, minZ, maxZ;
	float bX, bY, bZ;
	std::string path;
};