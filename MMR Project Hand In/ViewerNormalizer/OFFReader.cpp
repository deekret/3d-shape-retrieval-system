// basic file operations
#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>
#include <string>
#include <tuple>
#include "OFFReader.h"
#include "Grid.h"
#include "Renderer.h"
#include "FilterItem.h"
using namespace std;

typedef struct Vertex
{
	float x, y, z;					/* the usual 3-space position of a vertex */
} Vertex;

typedef struct Face
{
	unsigned char nverts;			/* number of vertex indices in list */
	vector<int> verts;              /* vertex index list */
} Face;


std::tuple<Grid*, FilterItem> openFile(string fileName)
{
	FilterItem fi;
	Grid* grid;

	int vertex_count, faces_count;

	ifstream first_file(fileName);

	if (first_file) {
		string line;

		/*while (getline(first_file, line)) {
			cout << line << endl;
		}*/

		if (first_file.good()) {
			if (getline(first_file, line));
			//cout << line << endl;
		}
		if (line == "ply") {										//If the file is in the .ply format, we use a differnet function
			
			fileName = readPlyFile(fileName);

		}
		first_file.close();

		ifstream offFile(fileName);
		getline(offFile, line);
		line = "";

		string counterLine;
		getline(offFile, counterLine);
		vector<string> string_vector;
		string_vector = split(counterLine, ' ');


		vertex_count = stoi(string_vector[0]);
		faces_count = stoi(string_vector[1]);

		fi.numVertices = stoi(string_vector[0]);
		fi.numFaces = stoi(string_vector[1]);
		if (fileName.find("/") <= size(fileName))
		{
			string folder = fileName.substr(fileName.find("/") + 1);
			fi.path = folder.substr(folder.find("/") + 1);
		}
		else
		{
			string folder = fileName.substr(fileName.find("\\") + 1);
			fi.path = folder.substr(folder.find("\\") + 1);
		}

		int i, j;
		grid = new Grid(fi.numVertices, fi.numFaces);

		float mix = FLT_MAX, miy = FLT_MAX, miz = FLT_MAX;
		float max = FLT_MIN, may = FLT_MIN, maz = FLT_MIN;
		float barix = 0, bariy = 0, bariz = 0;

		for (j = 0; j < fi.numVertices; j++)
		{
			getline(offFile, line);
			vector<string> vert;
			string_vector = split(line, ' ');
			Vertex v;
			vector<float> V;
			v.x = stof(string_vector[0]);
			v.y = stof(string_vector[1]);
			v.z = stof(string_vector[2]);
			V.push_back(v.x);
			V.push_back(v.y);
			V.push_back(v.z);
			grid->setPoint(j, V);

			max = std::max(v.x, max);
			may = std::max(v.y, may);
			maz = std::max(v.z, maz);
			
			mix = std::min(v.x, mix);
			miy = std::min(v.y, miy);
			miz = std::min(v.z, miz);

			barix += v.x;
			bariy += v.y;
			bariz += v.z;
		}

		fi.bX = barix / fi.numVertices;
		fi.bY = bariy / fi.numVertices;
		fi.bZ = bariz / fi.numVertices;

		fi.maxX = max;
		fi.maxY = may;
		fi.maxZ = maz;
		
		fi.minX = mix;
		fi.minY = miy;
		fi.minZ = miz;

		bool allTriangles = true, allQuads = true;

		for (j = 0; j < fi.numFaces; j++)
		{
			string line;
			getline(offFile, line);
			vector<string> vert;
			string_vector = split(line, ' ');
			Face face;

			int vertex_num = stoi(string_vector[0]);

			allTriangles &= vertex_num == 3;
			allQuads &= vertex_num == 4;

			for (int i = 0; i < vertex_num; i++) {
				face.verts.push_back(stoi(string_vector[i + 1]));
			}
			grid->setCell(j, face.verts);
		}

		if (allTriangles)
			fi.typeOfFace = "triangles";
		else if (allQuads)
			fi.typeOfFace = "quads";
		else
			fi.typeOfFace = "mix";

		return make_tuple(grid, fi);
	}


	else {
		cout << "FILE DOES NOT EXIST!" << endl;
		return make_tuple(grid, fi);
	}

}

FilterItem scanFile(string fileName)
{
	FilterItem fi;

	int vertex_count, faces_count;

	ifstream first_file(fileName);

	if (first_file) {
		string line;

		if (first_file.good()) {
			if (getline(first_file, line));
			//cout << line << endl;
		}
		if (line == "ply") {										//If the file is in the .ply format, we use a differnet function

			readPlyFile(fileName);

		}
		first_file.close();

		ifstream offFile(fileName);
		getline(offFile, line);
		line = "";

		string counterLine;
		getline(offFile, counterLine);
		vector<string> string_vector;
		string_vector = split(counterLine, ' ');

		fi.numVertices = stoi(string_vector[0]);
		fi.numFaces = stoi(string_vector[1]);

		if (fileName.find("/") <= size(fileName))
		{
			string folder = fileName.substr(fileName.find("/") + 1);
			fi.path = folder.substr(folder.find("/") + 1);
		}
		else
		{
			string folder = fileName.substr(fileName.find("\\") + 1);
			fi.path = folder.substr(folder.find("\\") + 1);
		}

		int i, j;

		float mix = FLT_MAX, miy = FLT_MAX, miz = FLT_MAX;
		float max = FLT_MIN, may = FLT_MIN, maz = FLT_MIN;
		float barix = 0, bariy = 0, bariz = 0;

		for (j = 0; j < fi.numVertices; j++)
		{
			getline(offFile, line);
			string_vector = split(line, ' ');
			Vertex v;
			v.x = stof(string_vector[0]);
			v.y = stof(string_vector[1]);
			v.z = stof(string_vector[2]);

			max = std::max(v.x, max);
			may = std::max(v.y, may);
			maz = std::max(v.z, maz);

			mix = std::min(v.x, mix);
			miy = std::min(v.y, miy);
			miz = std::min(v.z, miz);

			barix += v.x;
			bariy += v.y;
			bariz += v.z;
		}

		fi.bX = barix / fi.numVertices;
		fi.bY = bariy / fi.numVertices;
		fi.bZ = bariz / fi.numVertices;

		fi.maxX = max;
		fi.maxY = may;
		fi.maxZ = maz;

		fi.minX = mix;
		fi.minY = miy;
		fi.minZ = miz;

		bool allTriangles = true, allQuads = true;

		for (j = 0; j < fi.numFaces; j++)
		{
			string line;
			getline(offFile, line);
			int vertex_num = stoi(line.substr(0, 1));

			allTriangles &= vertex_num == 3;
			allQuads &= vertex_num == 4;
		}

		if (allTriangles)
			fi.typeOfFace = "triangles";
		else if (allQuads)
			fi.typeOfFace = "quads";
		else
			fi.typeOfFace = "mix";

		return fi;
	}


	else {
		cout << "FILE DOES NOT EXIST!" << endl;
		return fi;
	}
}

string readPlyFile(string fileName) {

	int vertex_count = 0;
	int faces_count = 0;

	ifstream ply_file(fileName);
	string new_name = fileName.substr(0, fileName.find(".ply"));
	ofstream off_file(new_name + ".off");

	string ply_line;

	while (getline(ply_file, ply_line)) {										//Go through the header line by line

		vector<string> string_vector;
		string_vector = split(ply_line, ' ');

		if (string_vector[0] == "end_header") break;
		else if (string_vector[0] == "element") {
			if (string_vector[1] == "vertex") {
				vertex_count = stoi(string_vector[2]);
			}
			else if (string_vector[1] == "face") {
				faces_count = stoi(string_vector[2]);
			}
		}
	
	}
	off_file << "OFF" << endl;
	off_file << to_string(vertex_count) + " " + to_string(faces_count)  + " 0" << endl;
	while (getline(ply_file, ply_line)) {
		off_file << ply_line << endl;
	}

	off_file.close();
	ply_file.close();

	return (new_name + ".off");

}

bool file_exists(const char* fileName) {
	ifstream infile(fileName);
	return infile.good();
}

vector<string> split(string str, char delimiter) {
	stringstream ss(str);
	vector<string> stringVec;
	string stringPiece;

	while (getline(ss, stringPiece, delimiter)) {
		stringVec.push_back(stringPiece);
	}
	return stringVec;
}

