#include "DungeonGenerator.h"

#include <algorithm>
#include <glad/glad.h>

#include <random>
#include <algorithm>
#include <numeric>
#include <ranges>

bool CompareX(Point& p1, Point& p2) {
	return p1.x < p2.x;
}

// Generates a random set of points on the screen
void DungeonGenerator::Init(unsigned size)
{
	std::srand(std::time(nullptr));

	m_DontUseEdge.resize(m_EdgeBoolVecSize);

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

	// Give each a tiny nudge randomly along the x or y dir
	for (auto& room : m_Rooms) {
		float xOffset = static_cast<float>(std::rand()) / static_cast<float>(RAND_MAX);
		float yOffset = static_cast<float>(std::rand()) / static_cast<float>(RAND_MAX);

		float xDir = static_cast<float>(std::rand()) / static_cast<float>(RAND_MAX);
		float yDir = static_cast<float>(std::rand()) / static_cast<float>(RAND_MAX);

		if (xDir < 0.5f) xOffset = -xOffset;
		if (yDir < 0.5f) yOffset = -yOffset;

		room.m_Position += glm::vec2(xOffset * 0.0001, yOffset * 0.0001);
	}

	// Run until all rooms are separated
	while (roomsNeedingSeparation) {

		for (auto& room : m_Rooms) {
			int neighborCount = 0;
			glm::vec2 v(0);

			// For each room
			for (auto& otherRoom : m_Rooms)
			{
				// If same room - skip
				if (&room == &otherRoom) continue;

				float roomPosDiff = glm::length(room.m_Position - otherRoom.m_Position);

				// If rooms are too close
				if (roomPosDiff < ((room.m_Width + otherRoom.m_Width) * 1.7f) ||
					roomPosDiff < ((room.m_Height + otherRoom.m_Height) * 1.7f))
				{
					neighborCount++;
					v += otherRoom.m_Position - room.m_Position;
				}
			}

			// If no separating vel - far away enough from every other room
			if (v == glm::vec2(0))
			{
				if (roomsNeedingSeparation > 0) roomsNeedingSeparation--;
				room.m_Velocity = glm::vec2(0);
				//continue;
			}
			else
			{
				v *= -1;
				v /= neighborCount;
				v = glm::normalize(v);

				room.m_Velocity = v * 10.0f;
				room.m_Position += (room.m_Velocity * 0.01f);
			}
		}
	}

	for (auto& room : m_Rooms) room.UpdateCorners();
}

Triangle DungeonGenerator::GenerateSuperTriangle()
{
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

	Point pA(x1, y1); pA.isSupra = true;
	Point pB(x2, y2); pB.isSupra = true;
	Point pC(x3, y3); pC.isSupra = true;

	// Arrange in anti-clockwise order
	float temp = (pB.x - pA.x) * (pC.y - pA.y) - (pC.x - pA.x) * (pB.y - pA.y);
	bool isCounterClockwise = temp > 0;
	Point vA = pA;
	Point vB = isCounterClockwise ? pB : pC;
	Point vC = isCounterClockwise ? pC : pB;

	// Place edges into the main vector. 
	// Edges will now be added in anti-clockwise order
	m_Edges.emplace_back(vA, vB);
	m_Edges.emplace_back(vB, vC);
	m_Edges.emplace_back(vC, vA);

	// Edges at indices 0, 1, and 2 of the
	// m_Edges vec make up this tri
	Triangle superTri(0, 1, 2);
	CalculateCircumCenter(superTri, vA, vB, vC);

	// Anti-clockwise order should be maintained
	return superTri;
}

// Bowyer-Watson's algorithm https://www.tinyurl.com/488kubkh
void DungeonGenerator::Triangulate()
{
	// Define a super/supra-triangle that surrounds all the points. 
	// Add them to the end of the list and mark them as isSupra (=true) 
	if (m_currentStepPointIndex == 0) {
		for (unsigned index = 0; index < m_Rooms.size(); index++) {
			m_Points.emplace_back(
				m_Rooms[index].m_Position.x,
				m_Rooms[index].m_Position.y,
				index);
		}
		// Sort the points in ascending order of x-coordinate
		std::sort(m_Points.begin(), m_Points.end(), CompareX);
		m_Triangles.push_back(GenerateSuperTriangle()); 

		size_t pointCount = m_Points.size();
		for (int index = 0; index < pointCount; index++) {
			m_Points[index].id = index;
		}
	}
	
	// Loop through points
	size_t pointsListSize = m_Points.size();
	for (unsigned pointIndex = m_currentStepPointIndex; 
		pointIndex < pointsListSize; pointIndex++) {

		Point& p = m_Points[pointIndex];

		// if is supra, we've reached the end. Leave
		if (p.isSupra) break;

		// indices of all bad edges
		std::vector<unsigned> badEdges;
		
		//	Examine the list of all triangles formed so far. 
		for (int triangleIndex = 0; triangleIndex < m_Triangles.size(); triangleIndex++) {

			auto& triangle = m_Triangles[triangleIndex];

			// For each triangle which is flagged as incomplete we must do distance checks
			if (triangle.status == Status::INCOMPLETE) {

				float R2 = triangle.circumRadiusSquared;

				float dx = (triangle.circumCenter.x - p.x);
				float dx2 = dx * dx;


				if (dx2 >= R2) {
					//	If Dx2 >= R2, then the circumclrcle for this triangle
					//	cannot be Intersected by any of the remaining points
					//	Flag this triangle as complete and do not execute
					//	steps 8 and 9
					triangle.status = Status::COMPLETE;
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
				if (badEdges[i] == badEdges[j]) {
					// We don't actually delete since it is potentially expensive
					// Just mark it as - Don't use.
					m_DontUseEdge[badEdges[i]] = true;
				}
			}
		}	
		 

		// Create new edges and triangles from those edges
		size_t badEdgeVecSize = badEdges.size();
		EdgeVec newEdges;
		
		// Loop though all indices of bad edges
		for (unsigned i = 0; i < badEdgeVecSize; i++) {
			int first, second;
			if (!m_DontUseEdge[badEdges[i]]) {

				Edge edge1(p, m_Edges[badEdges[i]].p1);
				auto edgeIter = std::find(newEdges.begin(), newEdges.end(), edge1);
				if (edgeIter == newEdges.end()) {
					newEdges.push_back(edge1); 
					first = m_Edges.size() + newEdges.size() - 1; 
				}
				else { 
					first = m_Edges.size() + std::distance(newEdges.begin(), edgeIter); }
				
				Edge edge2(p, m_Edges[badEdges[i]].p2);
				edgeIter = std::find(newEdges.begin(), newEdges.end(), edge2);
				if (edgeIter == newEdges.end()) {
					newEdges.push_back(edge2); 
					second = m_Edges.size() + newEdges.size() - 1; 
				}
				else { 
					second = m_Edges.size() + std::distance(newEdges.begin(), edgeIter); }

				m_Triangles.emplace_back(first, badEdges[i], second);

				Point pA = m_Edges[badEdges[i]].p1;
				Point pB = m_Edges[badEdges[i]].p2;
				Point pC = p;

				float temp = (pB.x - pA.x) * (pC.y - pA.y) - (pC.x - pA.x) * (pB.y - pA.y);
				bool isCounterClockwise = temp > 0;
				Point vA = pA;
				Point vB = isCounterClockwise ? pB : pC;
				Point vC = isCounterClockwise ? pC : pB;

				CalculateCircumCenter(m_Triangles.back(), vA, vB, vC);
			}
		}

		m_Edges.insert(m_Edges.end(), newEdges.begin(), newEdges.end());

		if (m_Edges.size() > m_DontUseEdge.size()) {
			m_EdgeBoolVecSize *= 2;
			m_DontUseEdge.resize(m_EdgeBoolVecSize);
		}

		if (m_StepWiseModeOn) {
			m_currentStepPointIndex = pointIndex + 1; break;
		}
	}

	if (m_StepWiseModeOn && m_currentStepPointIndex < m_Points.size()) {
		return;
	}

	//	Form the final triangulation by removing all triangles
	//	which have one or more of the supertriangle
	//	vertices.
	
	// TODO: Update this
	size_t triangleVecSize = m_Triangles.size();
	unsigned edgeAIndex, edgeBIndex, edgeCIndex;
	for (int index = 0; index < triangleVecSize; index++) {
		edgeAIndex = m_Triangles[index].edgeIndexA;
		edgeBIndex = m_Triangles[index].edgeIndexB;
		edgeCIndex = m_Triangles[index].edgeIndexC;
		if (m_Edges[edgeAIndex].p1.isSupra || m_Edges[edgeAIndex].p2.isSupra) {
			m_DontUseEdge[edgeAIndex] = true;
		}
		if (m_Edges[edgeBIndex].p1.isSupra || m_Edges[edgeBIndex].p2.isSupra) {
			m_DontUseEdge[edgeBIndex] = true;
		}
		if (m_Edges[edgeCIndex].p1.isSupra || m_Edges[edgeCIndex].p2.isSupra) {
			m_DontUseEdge[edgeCIndex] = true;
		}
	}
}


void DrawRoom(unsigned i, RoomSet& rooms) {
	auto room = rooms[i];
	glVertex2f(room.m_Corners[0].x, room.m_Corners[0].y);
	glVertex2f(room.m_Corners[1].x, room.m_Corners[1].y);

	glVertex2f(room.m_Corners[1].x, room.m_Corners[1].y);
	glVertex2f(room.m_Corners[2].x, room.m_Corners[2].y);

	glVertex2f(room.m_Corners[2].x, room.m_Corners[2].y);
	glVertex2f(room.m_Corners[3].x, room.m_Corners[3].y);

	glVertex2f(room.m_Corners[3].x, room.m_Corners[3].y);
	glVertex2f(room.m_Corners[0].x, room.m_Corners[0].y);
}

void DungeonGenerator::Display() {
	glLineWidth(1.0f);
	glPointSize(10.0f);
	glColor3f(1.0f, 1.0f, 1.0f);

	DisplayPoints();
	glBegin(GL_LINES);
	size_t edgeCount = m_FinalEdgesIndices.size();
	for (unsigned index = 0; index < edgeCount; index++) {
		auto& edge = m_Edges[m_FinalEdgesIndices[index]];

		DrawRoom(edge.p1.id, m_Rooms);
		DrawRoom(edge.p2.id, m_Rooms);

		glVertex2f(edge.p1.x, edge.p1.y);
		glVertex2f(edge.p2.x, edge.p2.y);
	}
	glEnd();
}

void DungeonGenerator::DisplayPoints() {

	glColor3f(1, 1, 1);

	glBegin(GL_POINTS);
	for (auto& p : m_Points) {
		glVertex2f(p.x, p.y);
	}
	glEnd();
}

void DungeonGenerator::CalculateCircumCenter(Triangle& tri, Point vA, Point vB, Point vC)
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


Point FindParent(Point point, 
	std::unordered_map<Point, Point>& parent) {
	if (parent[point] == point) {
		return point;
	}
	else {
		return FindParent(parent[parent[point]], parent);
	}

}

void DungeonGenerator::CreateMST()
{
	size_t edgeCount = m_Edges.size();
	for (unsigned index = 0; index < edgeCount; index++) {
		if (!m_DontUseEdge[index]) {
			m_ExtraEdgesIndices.push_back(index);
		}
	}

	auto pointCount = m_Points.size();

	std::vector<unsigned> mstEdges;

	// no. of edges = no. of nodes/vertices - 1
	mstEdges.resize(pointCount - 1);

	// Array to store parent of each point
	std::unordered_map<Point, Point> parents;
	
	// Initialize parents to self
	for (int index = 0; index < pointCount; index++) {
		parents[m_Points[index]] = m_Points[index];
	}

	unsigned count = 0, index = 0;
	while (count != pointCount - 1) {

		if (!m_DontUseEdge[index]) {
			Edge edge = m_Edges[index];

			// TODO: Find source parent and destination parent
			Point sourceParent = FindParent(edge.p1, parents);
			Point destParent = FindParent(edge.p2, parents);

			if (!(sourceParent == destParent)) {
				mstEdges[count++] = index;
				parents[sourceParent] = destParent;
			}
		}
		index++;
	}

	m_FinalEdgesIndices.insert(m_FinalEdgesIndices.end(),
		mstEdges.begin(), mstEdges.end());
}

void DungeonGenerator::AddBackExtraEdges()
{
	std::srand(std::time(nullptr));
	unsigned extraEdgeCount = m_ExtraEdgesIndices.size();
	for (int index = 0; index < extraEdgeCount; index++) {
		float randomNum = (std::rand() / static_cast<float>(RAND_MAX));
		if (randomNum < 0.3f) { 
			m_FinalEdgesIndices.push_back(m_ExtraEdgesIndices[index]);
		}
	}
}
