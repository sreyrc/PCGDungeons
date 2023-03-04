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

	///The float 0.866 is an arbitrary value determined for optimum supra triangle conditions.
	// TODO: Change back to 0.866f
	float x1 = xCen - 0.5f * dMax;
	float x2 = xCen + 0.5f * dMax;
	float x3 = xCen;

	float y1 = yCen - 0.5f * dMax;
	float y2 = yCen - 0.5f * dMax;
	float y3 = yCen + dMax;

	Point pointA(x1, y1); pointA.isSupra = true;
	Point pointB(x2, y2); pointB.isSupra = true;
	Point pointC(x3, y3); pointC.isSupra = true;

	m_Points.push_back(pointC);
	m_Points.push_back(pointA);
	m_Points.push_back(pointB);

	return Triangle(pointA, pointB, pointC);
}

// Bowyer-Watson's algorithm https://www.tinyurl.com/488kubkh
void DelaunayTriangulator::Triangulate()
{
	// Define a super/supra-triangle that surrounds all the points. 
	// Add them to the end of the list and mark them as isSupra (=true) 
	if (currStepPointIndex == 0) {
		// Sort the points in ascending order of x-coordinate
		std::sort(m_Points.begin(), m_Points.end(), CompareX);
		m_Triangles.push_back(GenerateSuperTriangle()); 
	}
	
	// Loop through points
	size_t pointsListSize = m_Points.size();
	for (unsigned pointIndex = currStepPointIndex; 
		pointIndex < pointsListSize; pointIndex++) {

		Point& p = m_Points[pointIndex];

		// if is supra, we've reached the end. Leave
		if (p.isSupra) break;

		EdgeVec edges;
		
		//	Examine the list of all triangles formed so far. 
		for (int j = 0; j < m_Triangles.size(); j++) {

			auto& triangle = m_Triangles[j];

			// For each triangle which is flagged as incomplete we must do distance checks
			if (triangle.status == Status::INCOMPLETE) {

				float R2 = triangle.GetCircumRadiusSquared();

				float dx = (triangle.GetCircumCenter().x - p.x);
				float dx2 = dx * dx;


				if (dx2 >= R2) {
					//	(7) If Dx2 >> -R~, then the circumclrcle for this triangle
					//	cannot be Intersected by any of the remaining points
					//	Flag this triangle as complete and do not execute
					//	steps 8 and 9
					triangle.status = Status::COMPLETE;
					continue;
				}

				//  Compute the square of the distance from the new
				//	point to the triangle clrcumcentre
				
				float dy = (triangle.GetCircumCenter().y - p.y);
				float d2 = dx2 + (dy * dy);
				 

				if (d2 < R2) {
					//	If D2 < R2, then the new point intersects the circum -
					//	circle for this triangle. Delete this triangle from the
					//	list of triangles and store the three pairs of
					//	vertices which define its edges on a list of edges. If
					//	D2 >= R2, then the new point hes on or outside the
					//	clrcumcircle for this triangle and the triangle remains
					//	unmodified

					edges.emplace_back(triangle.vA, triangle.vB);
					edges.emplace_back(triangle.vB, triangle.vC);
					edges.emplace_back(triangle.vC, triangle.vA);

					m_Triangles.erase(m_Triangles.begin() + j); j--;
				}
			}

		}	

		//	Loop over the list of edges and delete all edges which
		//	occur twice in the list. 
		std::vector<bool> dontUseEdge(edges.size());

		// TODO: Can we do better than a double loop?
		for (int i = 0; i < edges.size(); i++) {
			for (int j = i + 1; j < edges.size(); j++) {
				if (i == j) continue;
				if (edges[i] == edges[j]) {
					dontUseEdge[i] = true;
					dontUseEdge[j] = true;
				}
			}
		}
		 
		
		//	(11) Form the new triangles by matching the new point
		//	with each pair of vertices in the list of edges. The
		//	new point forms new triangles with each pair of
		//	vertices on the boundary of the polygon formed by
		//	the intersected triangles. Define each new triangle
		//	such that Its vertices are always listed in an anti -
		//	clockwise sequence and flag it as incomplete
		for (int i = 0; i < edges.size(); i++) {
			if(!dontUseEdge[i]) m_Triangles.emplace_back(Point(p.x, p.y), edges[i].p1, edges[i].p2);
		}

		if (m_StepWiseModeOn) {
			currStepPointIndex = pointIndex + 1; break;
		}
	}

	if (m_StepWiseModeOn && currStepPointIndex < m_Points.size() - 3) {
		return;
	}
	//	Form the final triangulatlon by removing all triangles
	//	which have one or more of the supertriangle
	//	vertices.
	
	for (int i = 0; i < m_Triangles.size(); i++) {
		if (m_Triangles[i].vA.isSupra ||
			m_Triangles[i].vB.isSupra ||
			m_Triangles[i].vC.isSupra) {

			m_Triangles.erase(m_Triangles.begin() + i); i--;
		}
	}
}

void DelaunayTriangulator::Display() {
	glLineWidth(1.5f);
	glPointSize(10.0f);
	glColor3f(1.0f, 1.0f, 1.0f);

	DisplayPoints();
	for (auto& tri : m_Triangles) {
		glBegin(GL_LINES); 
			glVertex2f(tri.vA.x, tri.vA.y);
			glVertex2f(tri.vB.x, tri.vB.y);

			glVertex2f(tri.vB.x, tri.vB.y);
			glVertex2f(tri.vC.x, tri.vC.y);

			glVertex2f(tri.vC.x, tri.vC.y);
			glVertex2f(tri.vA.x, tri.vA.y);
			
			//glColor3f(1, 1, 0);

			//glVertex2f(tri.vA.x, tri.vA.y);
			//glVertex2f(tri.GetCircumCenter().x, tri.GetCircumCenter().y);

			//glVertex2f(tri.vB.x, tri.vB.y);
			//glVertex2f(tri.GetCircumCenter().x, tri.GetCircumCenter().y);

			//glVertex2f(tri.vC.x, tri.vC.y);
			//glVertex2f(tri.GetCircumCenter().x, tri.GetCircumCenter().y);

			//glColor3f(1, 1, 1);
		glEnd();
	}
}

void DelaunayTriangulator::DisplayPoints() {

	glColor3f(1, 1, 1);

	glBegin(GL_POINTS);
	for (auto& p : m_Points) {
		glVertex2f(p.x, p.y);
	}
	glEnd();
}
