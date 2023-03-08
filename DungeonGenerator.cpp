// ---------------------------------------------------------------------------
// Written by Sreyash Raychaudhuri
// File Purpose: Implementation of DungeonGenerator class
// ---------------------------------------------------------------------------

#include <glad/glad.h>

#include <algorithm>
#include <ctime>

#include "DungeonGenerator.h"


/// @brief Point with smaller x coordinate is considered "smaller"
/// @param p1 First point
/// @param p2 Second point
/// @return true if p1.x < p2.x else false
bool CompareX(Point& p1, Point& p2) {
	return p1.x < p2.x;
}


void DungeonGenerator::Generate(unsigned numRooms)
{
	// Step 1: Randomly generate numRooms nummber of rooms
	Init(numRooms);

	// Step 2: Separate/Spread the rooms out
	Separate();

	// Step 3: Create a Delaunay Triangulation of the rooms (nodes)
	Triangulate();

	// Step 4: Create an MST from the output of the triangulation
	CreateMST();

	// Step 5: Add back some extra edges which were 
	// part of the triangulation but not part of the MST
	AddBackExtraEdges();
}


void DungeonGenerator::Init(unsigned size)
{
	std::srand(std::time(nullptr));

	m_EdgeStatus.resize(m_EdgeStatusVecSize);

	for (int index = 0; index < size; index++) {
		float roomWidthX = static_cast<float>(std::rand()) / static_cast<float>(RAND_MAX);
		float roomWidthY = static_cast<float>(std::rand()) / static_cast<float>(RAND_MAX);

		roomWidthX = std::clamp(roomWidthX, 0.4f, 0.8f);
		roomWidthY = std::clamp(roomWidthY, 0.4f, 0.8f);

		m_Rooms.emplace_back(roomWidthX * 100, roomWidthY * 100, index);
	}
}


void DungeonGenerator::Separate()
{
	unsigned roomsNeedingSeparation = m_Rooms.size();

	// All rooms are initially positioned at (0, 0)
	// To start separation there must be some difference in position between the rooms 
	// Therefore, we give each room a tiny nudge randomly along the x or y dir
	for (auto& room : m_Rooms) {
		float xOffset = static_cast<float>(std::rand()) / static_cast<float>(RAND_MAX);
		float yOffset = static_cast<float>(std::rand()) / static_cast<float>(RAND_MAX);

		float xDir = static_cast<float>(std::rand()) / static_cast<float>(RAND_MAX);
		float yDir = static_cast<float>(std::rand()) / static_cast<float>(RAND_MAX);

		if (xDir < 0.5f) xOffset = -xOffset;
		if (yDir < 0.5f) yOffset = -yOffset;

		room.m_Position += Point(xOffset * 0.0001, yOffset * 0.0001);
	}

	// Run until all rooms are separated
	while (roomsNeedingSeparation) {

		// For every room
		for (auto& room : m_Rooms) {
			int neighborCount = 0;
			Point v;

			// Check against other rooms
			for (auto& otherRoom : m_Rooms)
			{
				// If same room - skip
				if (&room == &otherRoom) continue;

				float roomPosDiff = room.m_Position.DistanceTo(otherRoom.m_Position);

				// If rooms are too close
				if (roomPosDiff < ((room.m_Width + otherRoom.m_Width) * 1.7f) ||
					roomPosDiff < ((room.m_Height + otherRoom.m_Height) * 1.7f))
				{
					neighborCount++;

					// Increase separating velocity
					v += (otherRoom.m_Position - room.m_Position);
				}
			}

			// If no separating vel - far away enough from every other room
			if (v == Point(0, 0))
			{
				if (roomsNeedingSeparation > 0) roomsNeedingSeparation--;
				room.m_Velocity = Point(0, 0);
			}
			else
			{
				// Adjust velocity
				v = v * -1;
				v = v / static_cast<float>(neighborCount);
				v.Normalize();

				room.m_Velocity = v * 10.0f;

				// Increment position depending on velocity
				room.m_Position += (room.m_Velocity * 0.01f);
			}
		}
	}

	// New positions have been set - now update the corners of the rooms
	for (auto& room : m_Rooms) room.UpdateCorners();
}


Triangle DungeonGenerator::GenerateSuperTriangle()
{
	// First find point bounds that surround all the points
	float minX = FLT_MAX;
	float minY = FLT_MAX;
	float maxX = -FLT_MAX;
	float maxY = -FLT_MAX;

	for (const auto p: m_Points)
	{
		if (minX > p.x) { minX = p.x; }
		if (minY > p.y) { minY = p.y; }
		if (maxX < p.x) { maxX = p.x; }
		if (maxY < p.y) { maxY = p.y; }
	}

	float dMax = std::max(maxX - minX, maxY - minY) * 3.0f;
	float xCen = (minX + maxX) * 0.5f;
	float yCen = (minY + maxY) * 0.5f;

	// The float 0.866 is an arbitrary value determined for optimum supra triangle conditions.
	float x1 = xCen - 0.866f * dMax;
	float x2 = xCen + 0.866f * dMax;
	float x3 = xCen;

	float y1 = yCen - 0.5f * dMax;
	float y2 = yCen - 0.5f * dMax;
	float y3 = yCen + dMax;

	// All these points are part of the supra triangle
	Point pA(x1, y1); pA.isSuper = true;
	Point pB(x2, y2); pB.isSuper = true;
	Point pC(x3, y3); pC.isSuper = true;

	// Place edges into the main vector. 
	// Edges will now be added in anti-clockwise order
	m_Edges.emplace_back(pA, pB);
	m_Edges.emplace_back(pB, pC);
	m_Edges.emplace_back(pC, pA);

	// Edges at indices 0, 1, and 2 of the
	// m_Edges vec make up this tri
	Triangle superTri(0, 1, 2);
	CalculateTriangleData(superTri, pA, pB, pC);

	// Anti-clockwise order should be maintained
	return superTri;
}


void DungeonGenerator::Triangulate()
{
	// Create a number of points based on the room positions
	for (unsigned index = 0; index < m_Rooms.size(); index++) {
		m_Points.emplace_back(
			m_Rooms[index].m_Position.x,
			m_Rooms[index].m_Position.y,
			index);
	}

	// Sort the points in ascending order of x-coordinate
	std::sort(m_Points.begin(), m_Points.end(), CompareX);

	// Define a super/supra-triangle that surrounds all the points. 
	// Add them to the end of the list and mark them as isSupra (=true) 
	m_Triangles.push_back(GenerateSuperTriangle()); 

	// Give each point and id.
	// id of point = it's index in the m_Points vec
	size_t pointCount = m_Points.size();
	for (int index = 0; index < pointCount; index++) {
		m_Points[index].m_ID = index;
	}
	
	// Loop through all points
	size_t pointsListSize = m_Points.size();
	for (unsigned pointIndex = 0; pointIndex < pointsListSize; pointIndex++) {

		Point& p = m_Points[pointIndex];

		// indices of all bad edges
		std::vector<unsigned> badEdges;
		
		//	Examine the list of all triangles formed so far. 
		for (int triangleIndex = 0; triangleIndex < m_Triangles.size(); triangleIndex++) {

			auto& triangle = m_Triangles[triangleIndex];

			// For each triangle which is flagged as incomplete we must do distance checks
			if (triangle.status == TriangleStatus::INCOMPLETE) {

				float R2 = triangle.circumRadiusSquared;

				float dx = (triangle.circumCenter.x - p.x);
				float dx2 = dx * dx;

				if (dx2 >= R2) {
					//	If Dx2 >= R2, then the circumclrcle for this triangle
					//	cannot be Intersected by any of the remaining points
					//	Flag this triangle as complete and do not execute
					//	steps 8 and 9
					triangle.status = TriangleStatus::COMPLETE;
					continue;
				}

				//  Compute the square of the distance from the new
				//	point to the triangle clrcumcentre
				
				float dy = (triangle.circumCenter.y - p.y);
				float d2 = dx2 + (dy * dy);
				 

				if (d2 < R2) {
					//	If D2 < R2, then the new point intersects the circum -
					//	circle for this triangle. Delete this triangle from the
					//	list of triangles and store the indices of edges (of the vector of edges) 
					//  of the triangle

					badEdges.emplace_back(triangle.edgeIndexA);
					badEdges.emplace_back(triangle.edgeIndexB);
					badEdges.emplace_back(triangle.edgeIndexC);

					// Remove this triangle
					m_Triangles.erase(m_Triangles.begin() + triangleIndex); 
					triangleIndex--;
				}
			}

		}		


		// Check for duplicate bad edge indices and remove them from the list
		for (unsigned i = 0; i < badEdges.size(); i++) {
			for (unsigned j = i + 1; j < badEdges.size(); j++) {
				if (i == j) continue;
				// If duplicate edge found
				if (badEdges[i] == badEdges[j]) {
					// Mark the edge as INVALID
					m_EdgeStatus[badEdges[i]] = EdgeStatus::INVALID;
				}
			}
		}	
		 

		// Create new edges and triangles from those edges
		size_t badEdgeCount = badEdges.size();
		EdgeVec newEdges;
		
		// Loop though all indices of bad edges
		for (unsigned i = 0; i < badEdgeCount; i++) {
			
			// two of the edge indices of the new triangle
			unsigned newTriEdge1, newTriEdge2;

			if (m_EdgeStatus[badEdges[i]] == EdgeStatus::AVAILABLE) {

				size_t edgeCount = m_Edges.size();
				
				// Create and edge from current point to one of the 
				// vertices of the bad edge
				Edge edge1(p, m_Edges[badEdges[i]].p1);

				// See if it's already been added
				auto edgeIter = std::find(newEdges.begin(), newEdges.end(), edge1);
				
				// If not, add it and set first index of new triangle accordingly
				if (edgeIter == newEdges.end()) {
					newEdges.push_back(edge1); 
					// Index of this edge = Size of m_Edges + the last index of newEdges
					// m_Edges will accommodate this edge when the new edges are appended to it
					newTriEdge1 = edgeCount + newEdges.size() - 1;
				}
				else { // else just set the index of first
					// Index = Size of m_Edges + index of newEdges where this edge was found
					newTriEdge1 = edgeCount + std::distance(newEdges.begin(), edgeIter); }
				
				// Do the same for the edge connecting current point and other point of the edge
				Edge edge2(p, m_Edges[badEdges[i]].p2);
				edgeIter = std::find(newEdges.begin(), newEdges.end(), edge2);
				if (edgeIter == newEdges.end()) {
					newEdges.push_back(edge2); 
					newTriEdge2 = edgeCount + newEdges.size() - 1;
				}
				else { 
					newTriEdge2 = edgeCount + std::distance(newEdges.begin(), edgeIter);
				}

				// Add this new triangle and calculate its data
				m_Triangles.emplace_back(newTriEdge1, badEdges[i], newTriEdge2);

				CalculateTriangleData(m_Triangles.back(), 
					m_Edges[badEdges[i]].p1, 
					m_Edges[badEdges[i]].p2,
					p);
			}
		}

		// Insert the set of new edges to the main edges vector
		m_Edges.insert(m_Edges.end(), newEdges.begin(), newEdges.end());

		// If EdgeStatus container too small, increase size
		if (m_Edges.size() > m_EdgeStatus.size()) {
			m_EdgeStatusVecSize *= 2;
			m_EdgeStatus.resize(m_EdgeStatusVecSize);
		}
	}


	//	Form the final triangulation by removing all triangles
	//	which have one or more of the super-triangle
	//	vertices.
	size_t triangleVecSize = m_Triangles.size();
	unsigned edgeAIndex, edgeBIndex, edgeCIndex;
	for (int index = 0; index < triangleVecSize; index++) {

		edgeAIndex = m_Triangles[index].edgeIndexA;
		edgeBIndex = m_Triangles[index].edgeIndexB;
		edgeCIndex = m_Triangles[index].edgeIndexC;
		
		// If either point of each edge of the triangle is part of supra triangle,
		// invlaidate this edge
		if (m_Edges[edgeAIndex].p1.isSuper || m_Edges[edgeAIndex].p2.isSuper) {
			m_EdgeStatus[edgeAIndex] = EdgeStatus::INVALID;
		}
		if (m_Edges[edgeBIndex].p1.isSuper || m_Edges[edgeBIndex].p2.isSuper) {
			m_EdgeStatus[edgeBIndex] = EdgeStatus::INVALID;
		}
		if (m_Edges[edgeCIndex].p1.isSuper || m_Edges[edgeCIndex].p2.isSuper) {
			m_EdgeStatus[edgeCIndex] = EdgeStatus::INVALID;
		}
	}
}


void DungeonGenerator::Display() {

	glClear(GL_COLOR_BUFFER_BIT);
	glClearColor(0, 0, 0, 1.0f);
	glLineWidth(1.5f);
	glPointSize(10.0f);
	glColor3f(1.0f, 1.0f, 0.0f);

	glBegin(GL_LINES);
	size_t edgeCount = m_Edges.size();

	// Display the dungeon edges
	for (unsigned index = 0; index < edgeCount; index++) {
		if (m_EdgeStatus[index] == EdgeStatus::IN_MST ||
			m_EdgeStatus[index] == EdgeStatus::EXTRA) {
			auto& edge = m_Edges[index];
			glVertex2f(edge.p1.x, edge.p1.y);
			glVertex2f(edge.p2.x, edge.p2.y);
		}
	}

	// Display the rooms
	glColor3f(1.0f, 1.0f, 1.0f);
	for (const auto& room : m_Rooms) {
		for (unsigned index = 0; index < 4; index++) {
			glVertex2f(room.m_Corners[index].x, room.m_Corners[index].y);
			glVertex2f(room.m_Corners[(index + 1) % 4].x, room.m_Corners[(index + 1) % 4].y);
		}
	}
	glEnd();

	// Display the room venter positions 
	// (initially used for debugging the graph algorithms)
	glBegin(GL_POINTS);
	for (const auto& room : m_Rooms) {
		glVertex2f(room.m_Position.x, room.m_Position.y);
	}
	glEnd();
}


void DungeonGenerator::CalculateTriangleData(Triangle& tri, Point vA, Point vB, Point vC)
{
	// Formulae: https://tinyurl.com/yvsw7s5f

	Point SqrA(powf(vA.x, 2), powf(vA.y, 2));
	Point SqrB(powf(vB.x, 2), powf(vB.y, 2));
	Point SqrC(powf(vC.x, 2), powf(vC.y, 2));

	float D = (vA.x * (vB.y - vC.y) + vB.x * (vC.y - vA.y) + vC.x * (vA.y - vB.y)) * 2.0f;
	float x = ((SqrA.x + SqrA.y) * (vB.y - vC.y) + (SqrB.x + SqrB.y) *
		(vC.y - vA.y) + (SqrC.x + SqrC.y) * (vA.y - vB.y)) / D;
	float y = ((SqrA.x + SqrA.y) * (vC.x - vB.x) + (SqrB.x + SqrB.y) *
		(vA.x - vC.x) + (SqrC.x + SqrC.y) * (vB.x - vA.x)) / D;

	tri.circumCenter = Point(x, y);

	// Square of Eucl. dist. between the CC and any of the vertices
	tri.circumRadiusSquared = ((vA.x - x) * (vA.x - x)) +
		((vA.y - y) * (vA.y - y));
}


/// @brief Find parent of given point at m_Points[pointIndex]
/// @param pointIndex Index of point in m_Points whose parent we wish to find
/// @param parent point parent vector
/// @return index of parent point
unsigned FindParent(unsigned pointIndex,
	std::vector<unsigned>& parent) {

	if (parent[pointIndex] == pointIndex) { 
		return pointIndex; 
	}
	else { 
		return FindParent(parent[pointIndex], parent); 
	}
}

void DungeonGenerator::CreateMST()
{
	size_t edgeCount = m_Edges.size();
	size_t pointCount = m_Points.size();

	std::vector<unsigned> mstEdges;

	// No. of edges = No. of nodes(vertices) - 1
	mstEdges.resize(pointCount - 1);

	// Array to store parent point index of each point
	// All indices of these points are to the m_Points vector
	std::vector<unsigned> parents(pointCount);
	
	// Initialize parents of points to self
	for (int index = 0; index < pointCount; index++) {
		parents[index] = index;
	}

	unsigned outputIndex = 0, inputIndex = 0;

	// Run the loop till we're at the end of the output vector
	while (outputIndex != pointCount - 1) {

		// Only if the edge is available for use
		if (m_EdgeStatus[inputIndex] == EdgeStatus::AVAILABLE) {
			Edge edge = m_Edges[inputIndex];

			// Find source parent and destination parent
			unsigned sourceParent = FindParent(edge.p1.m_ID, parents);
			unsigned destParent = FindParent(edge.p2.m_ID, parents);

			// If parents aren't the same - there's no cycle. Add this edge
			// to the MST and set parent source to dest.
			if (!(sourceParent == destParent)) {
				mstEdges[outputIndex++] = inputIndex;
				m_EdgeStatus[inputIndex] = EdgeStatus::IN_MST;
				parents[sourceParent] = destParent;
			}
		}

		// Go to next index in the m_Edges vec
		inputIndex++;
	}
}


void DungeonGenerator::AddBackExtraEdges()
{
	std::srand(std::time(nullptr));
	unsigned extraEdgeCount = m_Edges.size();
	for (int index = 0; index < extraEdgeCount; index++) {
		// Only if edge is availble for use
		if (m_EdgeStatus[index] == EdgeStatus::AVAILABLE) {
			// 30% chance that this edge may be marked as extra for use
			float randomNum = (std::rand() / static_cast<float>(RAND_MAX));
			if (randomNum < 0.3f) {
				m_EdgeStatus[index] = EdgeStatus::EXTRA;
			}
		}
	}
}