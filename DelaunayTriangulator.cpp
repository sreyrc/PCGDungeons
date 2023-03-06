#include "DelaunayTriangulator.h"

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
void DelaunayTriangulator::Init(unsigned size)
{
	std::vector<int> vX(size);
	std::vector<int> vY(size);

	std::random_device rd;
	std::iota(vX.begin(), vX.end(), 1);
	std::shuffle(vX.begin(), vX.end(), rd);

	std::iota(vY.begin(), vY.end(), 1);
	std::shuffle(vY.begin(), vY.end(), rd);

	for (unsigned index = 0; index < size; index++) {


		float newX = (((vX[index] - 1) * (2.0f)) / (float)(size - 1)) + -1.0f;
		float newY = (((vY[index] - 1) * (2.0f)) / (float)(size - 1)) + -1.0f;

		m_Points.emplace_back(newX, newY);
	}

	m_DontUseEdge.resize(m_EdgeBoolVecSize);
}

Triangle DelaunayTriangulator::GenerateSuperTriangle()
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
	// TODO: Change back to 0.866f
	float x1 = xCen - 0.5f * dMax;
	float x2 = xCen + 0.5f * dMax;
	float x3 = xCen;

	float y1 = yCen - 0.5f * dMax;
	float y2 = yCen - 0.5f * dMax;
	float y3 = yCen + dMax;

	Point pA(x1, y1); pA.isSupra = true;
	Point pB(x2, y2); pB.isSupra = true;
	Point pC(x3, y3); pC.isSupra = true;

	m_Points.push_back(pC);
	m_Points.push_back(pA);
	m_Points.push_back(pB);

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
void DelaunayTriangulator::Triangulate()
{
	// Define a super/supra-triangle that surrounds all the points. 
	// Add them to the end of the list and mark them as isSupra (=true) 
	if (m_currentStepPointIndex == 0) {
		// Sort the points in ascending order of x-coordinate
		std::sort(m_Points.begin(), m_Points.end(), CompareX);
		m_Triangles.push_back(GenerateSuperTriangle()); 
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

	if (m_StepWiseModeOn && m_currentStepPointIndex < m_Points.size() - 3) {
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

void DelaunayTriangulator::Display() {
	glLineWidth(1.0f);
	glPointSize(10.0f);
	glColor3f(1.0f, 1.0f, 1.0f);

	DisplayPoints();
	glBegin(GL_LINES);
	size_t edgeVecSize = m_Edges.size();
	for (unsigned i = 0; i < edgeVecSize; i++) {
		auto& edge = m_Edges[i];
		if (!m_DontUseEdge[i]) {
			glVertex2f(edge.p1.x, edge.p1.y);
			glVertex2f(edge.p2.x, edge.p2.y);
		}
	}
	glEnd();
}

void DelaunayTriangulator::DisplayPoints() {

	//glColor3f(1, 1, 1);

	//glBegin(GL_POINTS);
	//for (auto& p : m_Points) {
	//	glVertex2f(p.x, p.y);
	//}
	//glEnd();
}

void DelaunayTriangulator::CalculateCircumCenter(Triangle& tri, Point vA, Point vB, Point vC)
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
