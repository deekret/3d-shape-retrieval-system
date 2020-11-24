#include <GL/glew.h>
#include <GL/glut.h>
#include <GLFW/glfw3.h>
#include <filesystem>
#include <iostream>
#include <io.h>
#include <conio.h>
#include <filesystem>
#include <math.h>

#include <fstream>
#include <filesystem>
#include <tuple>
#include <stdio.h>
#include <string>
#include "OFFReader.h"
#include "Renderer.h"
#include "zpr.h"
using namespace std;
namespace fs = std::filesystem;
void loadFilter();
# define M_PI           3.14159265358979323846  /* pi */

string fileName;

int drawing_style = 0;
const int FILTER_SIZE = 250;


Grid* grid = 0;
Renderer renderer;

FilterItem* fis;

int index = 0;

float fov = 60;										                    //Perspective projection parameters
float z_near = -100;
float z_far = 100;

void draw()												                //Render the 3D mesh (GLUT callback function)
{
    renderer.draw(*grid);
}

void viewing(int W, int H)								//Window resize function, sets up viewing parameters (GLUT callback function)
{

    if (H == 0)
        H = 1;

    float ratio = W * 1.0 / W;

    glMatrixMode(GL_PROJECTION);						//1. Set the projection matrix
    glLoadIdentity();
    //gluPerspective(fov,float(W)/H,z_near,z_far);
    gluPerspective(45, ratio, 1, 100);
    glViewport(0, 0, W, H);
    glMatrixMode(GL_MODELVIEW);

    //glViewport(0,0,W,H);								//2. Set the viewport to the entire size of the rendering window
}

string getFileName(int index) {

    return ("m" + to_string(index));

}

void mouseclick(int button, int state, int x, int y)	//Callback for mouse click events. We use these to control the viewpoint interactively.
{
    int keys = glutGetModifiers();						//The mouse events are interpreted as follows:	
    if (keys & GLUT_ACTIVE_CTRL)							// 
        button = GLUT_MIDDLE_BUTTON;						// -left button + move:                         rotate viewpoint
    if (keys & GLUT_ACTIVE_SHIFT)							// -middle button (or left button+Control key): zoom in / out
        button = GLUT_RIGHT_BUTTON;							// -right buttom (or left button+Shift key):    translate viewpoint

    zprMouse(button, state, x, y);							//Use the ZPR library to manipulate the viewpoint. The library sets the modelview  
                                                          //OpenGL matrix based on the mouse events, thereby changing the current viewpoint.
}

void mousemotion(int x, int y)							//Callback for mouse move events. We use these to control the viewpoint interactively.
{
    zprMotion(x, y);										//Pass the current location of the mouse to the ZPR library, to change the viewpoint.

    glutPostRedisplay();									//After each mouse move event, ask GLUT to redraw our window, so we see the viewpoint change.
}

void mkdir(const char* dir)
{
    size_t size = strlen(dir) + 1;
    wchar_t* ndb = new wchar_t[size];
    size_t out;
    mbstowcs_s(&out, ndb, size, dir, size - 1);
    _wmkdir(ndb);
    delete[] ndb;
}

void outputFilter(string fileName)
{
    fstream filtout;
    filtout.open(fileName, ios::out);
    filtout << "sep=," << endl;
    filtout << "index,class,number of faces, number of vertices, type of faces, minimum X value, maximum X value, minimum Y value, maximum Y value, minimum Z value, maximum Z value, baricenter x coordinate, baricenter y coordinate, baricenter z coordinate, path, vertices verdict, volume" << endl;
    for (int i = 0; i < FILTER_SIZE; i++)
    {
        FilterItem fi = fis[i];

        if (fi.typeOfFace == "")
            continue;

        filtout << i;
        filtout << ",";
        filtout << fi.cls;
        filtout << ",";
        filtout << fi.numFaces;
        filtout << ",";
        filtout << fi.numVertices;
        filtout << ",";
        filtout << fi.typeOfFace;
        filtout << ",";
        filtout << fi.minX;
        filtout << ",";
        filtout << fi.maxX;
        filtout << ",";
        filtout << fi.minY;
        filtout << ",";
        filtout << fi.maxY;
        filtout << ",";
        filtout << fi.minZ;
        filtout << ",";
        filtout << fi.maxZ;
        filtout << ",";
        filtout << fi.bX;
        filtout << ",";
        filtout << fi.bY;
        filtout << ",";
        filtout << fi.bZ;
        filtout << ",";
        filtout << fi.path;
        filtout << ",";
        if (fi.numVertices < 1000)
        {
            filtout << "TOO FEW VERTICES";
            filtout << ",";
        }
        else if (fi.numVertices > 10000)
        {
            filtout << "TOO MANY VERTICES";
            filtout << ",";
        }
        else
        {
            filtout << ",";
        }
        filtout << (fi.maxX - fi.minX) * (fi.maxY - fi.minY) * (fi.maxZ - fi.minZ);
        filtout << endl;
    }

    cout << "Outputted!" << endl;
    filtout.close();
}

void scanFolder(string location)
{
    string folder;
    int i = 0;
    for (const auto& entry : fs::directory_iterator(location))
    {
        folder = entry.path().string();
        cout << folder << endl;
        for (const auto& fl : fs::directory_iterator(folder + "/"))
        {
            string file = fl.path().string();
            string suffix = ".off";
            if (!(file.size() >= suffix.size() && 0 == file.compare(file.size() - suffix.size(), suffix.size(), suffix)))
                continue;
            cout << file << endl;
            FilterItem fi = scanFile(file);
            int a = folder.find("/");
            if (a <= folder.size())
                fi.cls = folder.substr(folder.find("/") + 1);
            else
                fi.cls = folder.substr(folder.find("\\") + 1);
            fis[i] = fi;
            i++;
        }
    }
    cout << "Scanning complete!" << endl;
}


void keyboard(unsigned char c, int, int)					//Callback for keyboard events:
{
    switch (c)
    {
    case ' ':											    // space:   Toggle through the various drawing styles of the mesh renderer
    {
        index += 1;
        fileName = getFileName(index);
        cout << fileName << endl;//Grab the file, TODO: implement in a better way
        std::tuple<Grid*, FilterItem> tup = openFile("0/" + fileName + "/" + fileName + ".off");
        grid = std::get<0>(tup);
        fis[index] = std::get<1>(tup);

        break;
    }
    case 'D':
    case 'd':

        drawing_style = (++drawing_style) % 5;
        renderer.setDrawingStyle((Renderer::DRAW_STYLE)drawing_style);
        break;
    case 'o':
    {
        cout << "Please enter a name for the file" << endl;
        string fileName;
        cin >> fileName;
        outputFilter(fileName);
        break;
    }
    case 's':
    {
        cout << "Please enter a folder to scan" << endl;
        string fileName;
        cin >> fileName;
        scanFolder(fileName);
        break;
    }
    case 'l':
    {
        cout << "Loading from output file" << endl;
        loadFilter();
        break;
    }
    case 'n':
    {
        cout << "Normalizing current content..." << endl;

        string s = "Normalized_DB";
        mkdir(s.c_str());

        for (int i = 0; i < FILTER_SIZE; i++)
        {
            FilterItem fi = fis[i];

            if (fi.typeOfFace == "")
            {
                cout << i << endl;
                cout << "No known face type!" << endl;
                continue;
            }

            mkdir((s + "/" + fi.cls).c_str());
            tuple<Grid*, FilterItem> tup = openFile("Mechs/" + fi.cls + "/" + fi.path);

            Grid* g = get<0>(tup);
            float factor = max((fi.maxX - fi.minX), max((fi.maxY - fi.minY), (fi.maxZ - fi.minZ)));
            fstream fs;
            cout << s + "/" + fi.cls + "/" + fi.path << endl;
            fs.open(s + "/" + fi.cls + "/" + fi.path, ios::out);
            fs << "OFF" << endl;
            fs << fi.numVertices << " ";
            fs << fi.numFaces << " ";
            fs << 0 << endl;

            for (int j = 0; j < fi.numVertices; j++)
            {
                float* point = new float[3];
                g->getPoint(j, point);

                float fx = (point[0] - fi.bX) / factor;
                float fy = (point[1] - fi.bY) / factor;
                float fz = (point[2] - fi.bZ) / factor;

                fs << fx << " " << fy << " " << fz << endl;
            }

            //Write the new points into the file.
            for (int j = 0; j < g->numPoints(); j++)
            {
                float* point = new float[3];
                g->getPoint(j, point);

                float fx = point[0];
                float fy = point[1];
                float fz = point[2];

                //fs << fx << " " << fy << " " << fz << endl;
            }



            for (int j = 0; j < fi.numFaces; j++)
            {
                int* indices = new int[3];
                g->getCell(j, indices);

                fs << "3 " << indices[0] << " " << indices[1] << " " << indices[2] << endl;
            }

            fs.close();
            delete g;
        }
        cout << "Normalizing complete!" << endl;

        scanFolder("Normalized_DB");
        outputFilter("FilterOutput_after.csv");
        break;
    }
    case 't': {

        grid->PCARotation();
        grid->momentTest();
        
        float angle;
        float x, y, z;

        x = grid->eigenVec1[0];
        y = grid->eigenVec1[1];
        z = grid->eigenVec1[2];
        angle = acos(x / sqrt(pow(x, 2) + pow(y, 2) + pow(z, 2)));
        angle *= (180 / M_PI);
        cout << angle << endl;

        x = grid->eigenVec2[0];
        y = grid->eigenVec2[1];
        z = grid->eigenVec2[2];
        angle = acos(y / sqrt(pow(x, 2) + pow(y, 2) + pow(z, 2)));
        angle *= (180 / M_PI);
        cout << angle << endl;

        x = grid->eigenVec3[0];
        y = grid->eigenVec3[1];
        z = grid->eigenVec3[2];
        angle = acos(z / sqrt(pow(x, 2) + pow(y, 2) + pow(z, 2)));
        angle *= (180 / M_PI);
        cout << angle << endl;

        break;
    }
    case 'p': {

        scanFolder("Normalized_DB");

        fstream filtout;

        cout << "PCA Rotation..." << endl;

        string s = "Normalized_DB_PCA";
        mkdir(s.c_str());

        for (int i = 0; i < FILTER_SIZE; i++)
        {
            FilterItem fi = fis[i];

            if (fi.typeOfFace == "")
            {
                cout << i << endl;
                cout << "No known face type!" << endl;
                continue;
            }

            mkdir((s + "/" + fi.cls).c_str());
            tuple<Grid*, FilterItem> tup = openFile("Normalized_DB_Remeshed/" + fi.cls + "/" + fi.path);

            Grid* g = get<0>(tup);

            g->PCARotation();
            g->momentTest();
            grid->PCARotation();
            grid->momentTest();

            fstream fs;
            cout << s + "/" + fi.cls + "/" + fi.path << endl;
            fs.open(s + "/" + fi.cls + "/" + fi.path, ios::out);
            fs << "OFF" << endl;
            fs << fi.numVertices << " ";
            fs << fi.numFaces << " ";
            fs << 0 << endl;

            //Write the new points into the file.
            for (int j = 0; j < g->numPoints(); j++)
            {
                float* point = new float[3];
                g->getPoint(j, point);

                float fx = point[0];
                float fy = point[1];
                float fz = point[2];

                fs << fx << " " << fy << " " << fz << endl;
            }



            for (int j = 0; j < fi.numFaces; j++)
            {
                int* indices = new int[3];
                g->getCell(j, indices);

                fs << "3 " << indices[0] << " " << indices[1] << " " << indices[2] << endl;
            }

            fs.close();
            delete g;
        }
        cout << "PCA Rotation complete!" << endl;

        filtout.open("after_PCA.csv", ios::out);
        filtout << "sep=," << endl;
        filtout << "Index, angleX, angleY, angleZ, fX, fY, fZ" << endl;
        for (int i = 0; i < FILTER_SIZE; i++)
        {
            FilterItem fi = fis[i];

            if (fi.typeOfFace == "")
            {
                cout << i << endl;
                cout << "No known face type!" << endl;
                continue;
            }

            tuple<Grid*, FilterItem> tup = openFile("Normalized_DB_Remeshed_PCA/" + fi.cls + "/" + fi.path);
            Grid* g = get<0>(tup);
            g->computeEigenvectors();
            g->momentTest();

            if (fi.typeOfFace == "")
                continue;

            float angle;
            float x, y, z;

            filtout << i;
            filtout << ",";

            x = g->eigenVec1[0];
            y = g->eigenVec1[1];
            z = g->eigenVec1[2];
            angle = acos(x / sqrt(pow(x,2) + pow(y,2) + pow(z,2)));
            angle *= (180/ M_PI);
            filtout << angle;
            filtout << ",";

            x = g->eigenVec2[0];
            y = g->eigenVec2[1];
            z = g->eigenVec2[2];
            angle = acos(y / sqrt(pow(x, 2) + pow(y, 2) + pow(z, 2)));
            angle *= (180 / M_PI);
            filtout << angle;
            filtout << ",";

            x = g->eigenVec3[0];
            y = g->eigenVec3[1];
            z = g->eigenVec3[2];
            angle = acos(z / sqrt(pow(x, 2) + pow(y, 2) + pow(z, 2)));
            angle *= (180 / M_PI);
            filtout << angle;
            filtout << ",";


            filtout << g->f0;
            filtout << ",";
            filtout << g->f1;
            filtout << ",";
            filtout << g->f2;

            filtout << endl;
        }

        cout << "Outputted!" << endl;
        filtout.close();

        break;
    }
    }
    glutPostRedisplay();
}

void loadFilter()
{
    fstream filtin;
    string line;
    filtin.open("FilterOutput_after.csv", ios::in);
    if (filtin)
    {
        getline(filtin, line);
        getline(filtin, line);
        while (!filtin.eof())
        {
            getline(filtin, line);
            vector<string> vec = split(line, ',');
            if (size(vec) == 0)
                continue;
            FilterItem fi;
            fi.cls = vec[1];
            fi.numFaces = stoi(vec[2]);
            fi.numVertices = stoi(vec[3]);
            fi.typeOfFace = vec[4];
            fi.minX = stof(vec[5]);
            fi.maxX = stof(vec[6]);
            fi.minY = stof(vec[7]);
            fi.maxY = stof(vec[8]);
            fi.minZ = stof(vec[9]);
            fi.maxZ = stof(vec[10]);
            fi.bX = stof(vec[11]);
            fi.bY = stof(vec[12]);
            fi.bZ = stof(vec[13]);
            fi.path = vec[14];
            fis[stoi(vec[0])] = fi;
        }
        filtin.close();
    }
    else
    {
        std::cout << "No previous filter output found" << endl;
    }
}

void completeNormalization() {
    scanFolder("DB");

    cout << "Normalizing current content..." << endl;
    string s = "Normalized_DB";
    mkdir(s.c_str());

    for (int i = 0; i < FILTER_SIZE; i++)
    {
        FilterItem fi = fis[i];

        if (fi.typeOfFace == "")
        {
            cout << i << endl;
            cout << "No known face type!" << endl;
            continue;
        }

        mkdir((s + "/" + fi.cls).c_str());
        tuple<Grid*, FilterItem> tup = openFile("DB/" + fi.cls + "/" + fi.path);

        Grid* g = get<0>(tup);
        float factor = max((fi.maxX - fi.minX), max((fi.maxY - fi.minY), (fi.maxZ - fi.minZ)));
        fstream fs;
        cout << s + "/" + fi.cls + "/" + fi.path << endl;
        fs.open(s + "/" + fi.cls + "/" + fi.path, ios::out);
        fs << "OFF" << endl;
        fs << fi.numVertices << " ";
        fs << fi.numFaces << " ";
        fs << 0 << endl;

        for (int j = 0; j < fi.numVertices; j++)
        {
            float* point = new float[3];
            g->getPoint(j, point);

            float fx = (point[0] - fi.bX) / factor;
            float fy = (point[1] - fi.bY) / factor;
            float fz = (point[2] - fi.bZ) / factor;

            //Write the new points into the file.
            fs << fx << " " << fy << " " << fz << endl;
        }

        //Write the new points into the file.
        for (int j = 0; j < g->numPoints(); j++)
        {
            float* point = new float[3];
            g->getPoint(j, point);

            float fx = point[0];
            float fy = point[1];
            float fz = point[2];

            //fs << fx << " " << fy << " " << fz << endl;
        }

        for (int j = 0; j < fi.numFaces; j++)
        {
            int* indices = new int[3];
            g->getCell(j, indices);
            fs << "3 " << indices[0] << " " << indices[1] << " " << indices[2] << endl;
        }
        fs.close();
        delete g;
    }

    cout << "Normalizing phase 1 complete!" << endl;
    scanFolder("Normalized_DB");
    //outputFilter("FilterOutput_after.csv");
    //scanFolder("Normalized_DB");
    fstream filtout;
    cout << "PCA Rotation..." << endl;
    s = "Normalized_DB";
    mkdir(s.c_str());
    for (int i = 0; i < FILTER_SIZE; i++)
    {
        FilterItem fi = fis[i];

        if (fi.typeOfFace == "")
        {
            cout << i << endl;
            cout << "No known face type!" << endl;
            continue;
        }

        mkdir((s + "/" + fi.cls).c_str());
        tuple<Grid*, FilterItem> tup = openFile("Normalized_DB/" + fi.cls + "/" + fi.path);

        Grid* g = get<0>(tup);

        g->PCARotation();
        g->momentTest();
        //grid->PCARotation();
        //grid->momentTest();

        fstream fs;
        cout << s + "/" + fi.cls + "/" + fi.path << endl;
        fs.open(s + "/" + fi.cls + "/" + fi.path, ios::out);
        fs << "OFF" << endl;
        fs << fi.numVertices << " ";
        fs << fi.numFaces << " ";
        fs << 0 << endl;

        //Write the new points into the file.
        for (int j = 0; j < g->numPoints(); j++)
        {
            float* point = new float[3];
            g->getPoint(j, point);

            float fx = point[0];
            float fy = point[1];
            float fz = point[2];

            fs << fx << " " << fy << " " << fz << endl;
        }

        for (int j = 0; j < fi.numFaces; j++)
        {
            int* indices = new int[3];
            g->getCell(j, indices);

            fs << "3 " << indices[0] << " " << indices[1] << " " << indices[2] << endl;
        }

        fs.close();
        delete g;
    }
    cout << "PCA Rotation complete!" << endl;
}

int main(int argc, char* argv[])
{

    cout << "Press 'v' to view a off or ply file." << endl;
    cout << "Press 'n' to do a full normalization on a data base." << endl;
    char  input = _getch();

    fis = new FilterItem[FILTER_SIZE];
    loadFilter();

    if (input == 'v') {

        string input;
        cout << "Please specify the file you want to view:" << endl;
        cin >> input;

        std::tuple<Grid*, FilterItem> tup = openFile(input);
        grid = std::get<0>(tup);
        fis[index] = std::get<1>(tup);

        grid->computeFaceNormals();
        grid->computeVertexNormals();

        glutInit(&argc, argv);								                //Initialize the GLUT toolkit
        glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE | GLUT_DEPTH);
        //Ask GLUT to create next windows with a RGB framebuffer and a Z-buffer too
        glutInitWindowSize(500, 500);							            //Tell GLUT how large are the windows we want to create next
        glutCreateWindow(fileName.c_str());	                                //Create our window

        glutMouseFunc(mouseclick);							                //Bind the mouse click and mouse drag (click-and-move) events to callbacks. This allows us
        glutMotionFunc(mousemotion);
        glutKeyboardFunc(keyboard);
        glutDisplayFunc(draw);
        //glutReshapeFunc(viewing);
        glutMainLoop();
    }
    if (input == 'n') {
        completeNormalization();
    }
    return 0;
}