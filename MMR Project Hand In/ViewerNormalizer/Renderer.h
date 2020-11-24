#pragma once

#include "Renderer.h"

class Grid;



class Renderer
{
public:

	enum DRAW_STYLE
	{
		DRAW_GRID = 0,						//Draw the grid only. See the SimpleRenderer class.
		DRAW_POINTS,						//Draw points only
		DRAW_C0_CELLS,						//Draw cells, using flat shading
		DRAW_C1_CELLS,						//Draw cells, using smooth shading
		DRAW_GRID_CELLS						//Draw cells with grid ontop
	};

	Renderer() : draw_style(DRAW_GRID) {}

	void			draw(Grid&);

	void			setDrawingStyle(DRAW_STYLE s)
	{
		draw_style = s;
	}
	void drawAxis();

protected:

	void drawGrid(Grid&);
	void drawPoints(Grid&);
	void drawC0Cells(Grid&);
	void drawC1Cells(Grid&);

	DRAW_STYLE		draw_style;
};





#pragma once
