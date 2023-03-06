#pragma once

#include <vector>

#include "Primitives.h"

typedef std::vector<Point> PointVec;
typedef std::vector<Edge> EdgeVec;
typedef std::vector<Edge*> EdgePtrVec;
typedef std::vector<Triangle> TriangleVec;

typedef std::vector<bool> BoolVec;


class DelaunayTriangulator
{
public:
	void Init(unsigned size);
	Triangle GenerateSuperTriangle();
	void Triangulate();
	void Display();
	void DisplayPoints();
	void CalculateCircumCenter(Triangle& tri, 
		Point vA, Point vB, Point vC);
	bool StepWiseModeOn() { return m_StepWiseModeOn; }

private:
	PointVec m_Points;
	TriangleVec m_Triangles;
	bool m_StepWiseModeOn = false;
	unsigned m_currentStepPointIndex = 0;

	EdgeVec m_Edges;
	std::vector<bool> m_DontUseEdge;
	unsigned m_EdgeBoolVecSize = 100;
};

