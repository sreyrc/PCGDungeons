// ---------------------------------------------------------------------------
// Written by Sreyash Raychaudhuri
// File Purpose: Declaration of Room struct and DungeonGenerator class 
// DungeonGenerator procedurally generates a dungeon layout given a number of rooms 
// ---------------------------------------------------------------------------

#pragma once

#include <vector>
#include "Primitives.h"

/// @brief User-defined class that defines a room in the dungeon. 
/// Each room has a width, height and position, and corners.
/// Velocity is used for room separation
struct Room {

	/// @brief Width and height of the room
	float m_Width, m_Height;

	/// @brief Position and velocity of the room
	Point m_Position, m_Velocity;

	/// @brief Unique id defining the rooms
	unsigned m_ID;

	Room(float widthUnits, float heightUnits, int id)
		: m_Width (widthUnits * (1.0f / 360.0f)),
		m_Height (heightUnits * (1.0f / 360.0f)), m_ID(id) {

		m_Corners.resize(4);
		UpdateCorners();
	}

	/// @brief Update corners of the room based on current position
	void UpdateCorners() {
		m_Corners[0] = m_Position + Point(m_Height, m_Width);
		m_Corners[1] = m_Position + Point(m_Height, -m_Width);
		m_Corners[2] = m_Position + Point(-m_Height, -m_Width);
		m_Corners[3] = m_Position + Point(-m_Height, m_Width);
	}

	/// @brief Corner positions of the room
	std::vector<Point> m_Corners;
};

typedef std::vector<Point> PointVec;
typedef std::vector<Edge> EdgeVec;
typedef std::vector<Triangle> TriangleVec;
typedef std::vector<Room> RoomVec;


/// @brief A class that generates a dungeon 
/// layout given a number of rooms
class DungeonGenerator
{
public:
	/// @brief Generate a dungeon with the desired number of rooms  
	void Generate(unsigned numRooms);

	/// @brief Display the dungeon
	void Display();

private:
	/// @brief Generates rooms with random dimensions (with some clamping)
	/// @param size Number of rooms to be generated
	void Init(unsigned size);


	/// @brief Separate out the rooms using 
	/// a separation flocking algorithm: https://tinyurl.com/5h76ee97
	void Separate();


	/// @brief Generate a super triangle that surrounds 
	/// all the points to be triangulated.
	/// @return Super triangle surrounding all points
	Triangle GenerateSuperTriangle();


	/// @brief Delaunay Triangulation using
	/// Bowyer-Watson's algorithm: https://www.tinyurl.com/488kubkh.
	/// The algorithm in the paper returns the set of triangles in the triangulation -
	/// which means a lot of duplicate edges would be stored.
	/// Modified it such that only all the unique edges 
	/// that make up the triangulation are stored
	void Triangulate();

	/// @brief Calculate circumcircle and square of the circumradius of the triangle
	/// @param tri Triangle whose circumcenter is to be calculated
	/// @param vA First vertex of triangle tri
	/// @param vB Second vertex of triangle tri
	/// @param vC Third vertex of triangle tri
	void CalculateTriangleData(Triangle& tri,
		Point vA, Point vB, Point vC);
	

	/// @brief Find minimum spanning tree of the tirangulation output 
	/// edge-set using Kruskal's algorithm
	void CreateMST();


	/// @brief Randomly mark some of the AVAILABLE edges 
	/// as EXTRA for use of the dungeons
	void AddBackExtraEdges();

private:
	/// @brief All the rooms of the dungeon
	RoomVec m_Rooms;

	/// @brief List of all the points in the graph
	PointVec m_Points;

	/// @brief All the triangles for Delaunay Triangulation
	TriangleVec m_Triangles;

	/// @brief A list of all edges accumulated after running
	/// different algorithms
	EdgeVec m_Edges;

	/// @brief The statuses of every edge
	std::vector<EdgeStatus> m_EdgeStatus;

	/// @brief Size of the edge status vector
	unsigned m_EdgeStatusVecSize = 100;
};

