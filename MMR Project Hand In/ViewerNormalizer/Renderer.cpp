#include "Renderer.h"
#include "Grid.h"
#include <GL/glut.h>										//GLUT library

void Renderer::drawAxis() {
	glBegin(GL_LINES);

	glColor3f(1.0, 0.0, 0.0);
	glVertex3f(0.0, 0.0, 0.0);
	glVertex3f(4.0, 0.0, 0.0);

	glColor3f(0.0, 1.0, 0.0);
	glVertex3f(0.0, 0.0, 0.0);
	glVertex3f(0.0, 4.0, 0.0);

	glColor3f(0.0, 0.0, 1.0);
	glVertex3f(0.0, 0.0, 0.0);
	glVertex3f(0.0, 0.0, 4.0);
	glEnd();
}

void Renderer::draw(Grid& g_)								//Draw an unstructured grid, in the style indicated by 'drawing_style'
{
	glClearColor(1, 1, 1, 1);								//Clear the frame and depth buffers
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glEnable(GL_DEPTH_TEST);
	glEnable(GL_NORMALIZE);									//Ask OpenGL to automatically normalize normal-vectors when we scale the model
															//This is required to get correct shading for various model scales.

	Grid& g = (Grid&)g_;
	vector<Point3d> centroids;

	switch (draw_style)										//Render the grid, based on the current drawing style
	{														//Currently only implemented the grid version
	case DRAW_GRID:
		drawAxis();

		// //
		/*glColor3f(0.3, 0.3, 0.3);
		glEnable(GL_LIGHTING);
		glEnable(GL_LIGHT0);
		glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_TRUE);
		glShadeModel(GL_SMOOTH);
		glPointSize(3);
		centroids = g.getCellCentroids();
		glBegin(GL_POINTS);
		for (int i = 0; i < centroids.size(); i++) {
			glVertex3f(centroids[i].x, centroids[i].y, centroids[i].z);
		}*/
		// //

		glColor3f(1, 1, 1);
		glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
		drawGrid(g);
		glColor3f(0, 0, 0);
		glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
		drawGrid(g);
		break;
	case DRAW_POINTS:

		drawAxis();
		drawPoints(g);
		break;
	case DRAW_C0_CELLS:
		drawAxis();
		drawC0Cells(g);
		break;
	case DRAW_C1_CELLS:
		drawAxis();
		drawC1Cells(g);
		break;

	case DRAW_GRID_CELLS:

		drawAxis();
		glEnable(GL_POLYGON_OFFSET_FILL);
		glPolygonOffset(1.0, 2);

		drawC1Cells(g);
		glColor3f(1, 1, 1);
		glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
		drawGrid(g);
		glColor3f(0, 0, 0);
		glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
		drawGrid(g);
		break;
	}
	glutSwapBuffers();
	
}



void Renderer::drawGrid(Grid& g)							//Draw 'g' without shading
{

	glDisable(GL_LIGHTING);

	for (int i = 0; i < g.numCells(); ++i)					//Draw all cells in the grid
	{
		int   cell[10];										//

		int nvertices = g.getCell(i, cell);					//Get the 'nvertices' vertex-IDs of i-th cell
		if (nvertices != 3 && nvertices != 4)				//We only handle here drawing of triangle and quad cells.
			continue;										//It is quite simple to extend this to other cell types.

		if (nvertices == 3)									//Triangle cells:
			glBegin(GL_TRIANGLES);
		else //nvertices==4	
			glBegin(GL_QUADS);								//Quad cells:		

		for (int j = 0; j < nvertices; ++j)					//Render current cell
		{
			float x[3];
			g.getPoint(cell[j], x);							//Get vertex coordinates of j-th vertex of i-th cell		
			glVertex3f(x[0], x[1], x[2]);					//Pass this coordinate to OpenGL
		}
		glEnd();											//Finish rendering current cell
	}
}

void Renderer::drawPoints(Grid& g)	//Draw 'g' as a point cloud
{
	glColor3f(0.3, 0.3, 0.3);								//0. Base color for the grid points
	glEnable(GL_LIGHTING);
	glEnable(GL_LIGHT0);
	glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_TRUE);
	glShadeModel(GL_SMOOTH);
	glPointSize(2);										//1.   Draw all vertices in the grid as large points

	VectorAttributes& point_normals = g.getPointNormals();

	glBegin(GL_POINTS);									//2.   Draw vertices as an OpenGL point-set	

	for (int i = 0; i < g.numPoints(); ++i)
	{
		float v[3];
		g.getPoint(i, v);								//3. Get coordinates and normal of i-th vertex			
		float* point_normal = point_normals.getC0Vector(i);

		glNormal3f(point_normal[0], point_normal[1], point_normal[2]);
		glVertex3f(v[0], v[1], v[2]);						//2.2.  Pass normal, then coordinate, to OpenGL
	}
	glEnd();											//3.    Finish rendering the entire point-set
}

void Renderer::drawC0Cells(Grid& g)	//Draw 'g' with flat shading (one color per cell)
{
	glColor3f(0.3, 0.3, 0.3);								//2. Base color for the grid cells
	glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);			//3. Render cells filled
	glEnable(GL_LIGHTING);
	glEnable(GL_LIGHT0);
	glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_TRUE);
	glShadeModel(GL_FLAT);

	VectorAttributes& face_normals = g.getFaceNormals();

	for (int i = 0; i < g.numCells(); ++i)						//1. Draw all cells in the grid
	{
		int   cell[10];									//

		int nvertices = g.getCell(i, cell);				//2.   Get the 'nvertices' vertex-IDs of i-th cell
		if (nvertices != 3 && nvertices != 4)				//     We only handle here drawing of triangle and quad cells.
			continue;									//     It is quite simple to extend this to other cell types.

		if (nvertices == 3)								//     Triangle cells:
			glBegin(GL_TRIANGLES);
		else //nvertices==4	
			glBegin(GL_QUADS);							//     Quad cells:


		float* face_normal = face_normals.getC0Vector(i);//3.  Get normal of current cell, and pass it to OpenGL
		glNormal3f(face_normal[0], face_normal[1], face_normal[2]);

		for (int j = 0; j < nvertices; ++j)					//4.   Render current cell
		{
			float x[3];
			g.getPoint(cell[j], x);						//     Get vertex coordinates of j-th vertex of i-th cell		
			glVertex3f(x[0], x[1], x[2]);					//     Pass this coordinate to OpenGL
		}
		glEnd();                                        //5.   Finish rendering current cell
	}
}

void Renderer::drawC1Cells(Grid& g)	//Draw 'g' with smooth shading (bilinear interpolation of shading at vertices)
{
	glColor3f(0.3, 0.3, 0.3);								//2. Base color for the grid cells
	glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);			//3. Render cells filled
	glEnable(GL_LIGHTING);
	glEnable(GL_LIGHT0);
	glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_TRUE);
	glShadeModel(GL_SMOOTH);

	VectorAttributes& point_normals = g.getPointNormals();

	for (int i = 0; i < g.numCells(); ++i)						//1. Draw all cells in the grid 
	{
		int   cell[10];									//

		int nvertices = g.getCell(i, cell);				//2.   Get the 'nvertices' vertex-IDs of i-th cell
		if (nvertices != 3 && nvertices != 4)				//     We only handle here drawing of triangle and quad cells.
			continue;									//     It is quite simple to extend this to other cell types.

		if (nvertices == 3)								//     Triangle cells:
			glBegin(GL_TRIANGLES);
		else //nvertices==4	
			glBegin(GL_QUADS);							//     Quad cells:


		for (int j = 0; j < nvertices; ++j)					//3.   Render current cell
		{
			float x[3];									//     Get vertex coordinates of j-th vertex of i-th cell
			g.getPoint(cell[j], x);						//     Next, get vertex normal of the same vertex			
			float* point_normal = point_normals.getC0Vector(cell[j]);

			glNormal3f(point_normal[0], point_normal[1], point_normal[2]);
			glVertex3f(x[0], x[1], x[2]);					//     Pass the vertex normal, then the vertex coordinate, to OpenGL
		}
		glEnd();                                        //4.   Finish rendering current cell
	}
}








