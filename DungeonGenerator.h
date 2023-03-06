#pragma once

#include <vector>
#include <unordered_map>
#include <glm/glm.hpp>

#include "Primitives.h"

struct Room {

	float m_Width;
	float m_Height;

	glm::vec2 m_Position = glm::vec2(0);
	glm::vec2 m_Velocity = glm::vec2(0);

	unsigned ID;

	Room(float widthUnits, float heightUnits, int id)
		: m_Width (widthUnits * (1.0f / 360.0f)),
		m_Height (heightUnits * (1.0f / 360.0f)), ID(id) {

		m_Corners.resize(4);
		UpdateCorners();
	}

	void UpdateCorners() {
		m_Corners[0] = m_Position + glm::vec2(m_Height, m_Width);
		m_Corners[1] = m_Position + glm::vec2(m_Height, -m_Width);
		m_Corners[2] = m_Position + glm::vec2(-m_Height, -m_Width);
		m_Corners[3] = m_Position + glm::vec2(-m_Height, m_Width);
	}

	std::vector<glm::vec2> m_Corners;
};


typedef std::vector<Point> PointVec;
typedef std::vector<Edge> EdgeVec;
typedef std::vector<Triangle> TriangleVec;
typedef std::vector<Room> RoomSet;


class DungeonGenerator
{
public:
	void Init(unsigned size);
	void Separate();

	Triangle GenerateSuperTriangle();
	void Triangulate();
	void Display();
	void DisplayPoints();
	void CalculateCircumCenter(Triangle& tri, 
		Point vA, Point vB, Point vC);
	bool StepWiseModeOn() { return m_StepWiseModeOn; }

	void CreateMST();
	void AddBackExtraEdges();

private:
	RoomSet m_Rooms;

	PointVec m_Points;
	TriangleVec m_Triangles;
	bool m_StepWiseModeOn = false;
	unsigned m_currentStepPointIndex = 0;

	EdgeVec m_Edges;
	std::vector<unsigned> m_FinalEdgesIndices;
	std::vector<unsigned> m_ExtraEdgesIndices;
	std::vector<bool> m_DontUseEdge;
	unsigned m_EdgeBoolVecSize = 100;
};

