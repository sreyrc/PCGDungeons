#pragma once

#include <vector>

#include "Primitives.h"

typedef std::vector<Point> PointVec;
typedef std::vector<Edge> EdgeVec;
typedef std::vector<Triangle> TriangleVec;


class DelaunayTriangulator
{
public:
	void Init(unsigned size);
	Triangle GenerateSuperTriangle();
	void Triangulate();
	void Display();
	void DisplayPoints();

private:
	// TODO: Make var naming consistent
	PointVec m_Points;// m_SuperTriPoints;
	TriangleVec m_Triangles;
	bool m_StepWiseModeOn = true;
	unsigned currStepPointIndex = 0;
};

