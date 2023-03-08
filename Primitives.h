// ---------------------------------------------------------------------------
// Written by Sreyash Raychaudhuri
// File Purpose:  Definition of geometric primitives - Point, Edge and Triangle
// ---------------------------------------------------------------------------

#pragma once

/// @brief Point class - Acts like a Vector [X, Y]
class Point {
public:
	Point() : x(0.0f), y(0.0f), m_ID(INT_MAX) {}

	Point(float _x, float _y, int _id = INT_MAX, bool _isSupra = false)
		: x(_x), y(_y), m_ID(_id), isSuper(_isSupra) {}

	Point(const Point& p) {
		x = p.x;
		y = p.y;
		m_ID = p.m_ID;
		isSuper = p.isSuper;
	}

	float DistanceTo(const Point& p) const {
		return sqrtf(((x - p.x) * (x - p.x)) +
			((y - p.y) * (y - p.y)));
	}

	void Normalize() {
		float denom = sqrtf((x * x) + (y * y));
		x /= denom;
		y /= denom;
	}


	Point& operator=(const Point& p) {
		x = p.x;
		y = p.y;
		m_ID = p.m_ID;
		isSuper = p.isSuper;
		return (*this);
	}

	Point operator+(const Point& p) {
		Point pRes(x, y);
		pRes.x = x + p.x;
		pRes.y = y + p.y;
		return pRes;
	}

	Point operator-(const Point& p) {
		Point pRes(x, y);
		pRes.x = x - p.x;
		pRes.y = y - p.y;
		return pRes;
	}

	Point operator*(float factor) {
		Point pRes(x, y);
		pRes.x = x * factor;
		pRes.y = y * factor;
		return pRes;
	}

	Point operator/(float factor) {
		Point pRes(x, y);
		pRes.x = x / factor;
		pRes.y = y / factor;
		return pRes;
	}

	void operator+=(const Point& p) {
		x += p.x;
		y += p.y;
	}

	bool operator<(const Point& p) {
		return (x == p.x) ? ((y < p.y) ? true : false) :
			(x < p.x) ? true : false;
	}

	bool operator==(const Point& p) const {
		return (fabs(x - p.x) < FLT_EPSILON)
			&& (fabs(y - p.y) < FLT_EPSILON);
	}

	float x, y;
	unsigned m_ID;
	bool isSuper = false;
};


/// @brief Status of the edge:
/// AVAILABLE: Avaialable for use in algorithms
/// INVALID: Edge is not part of the triangulation
/// IN_MST: Edge is part of the MST
/// EXTRA: Edges not part of the MST but used in the dungeon
enum class EdgeStatus {
	AVAILABLE,
	INVALID,
	IN_MST,
	EXTRA
};

/// @brief Edge connects two Points
class Edge {
public:
	Edge() {}
	Edge(Point _p1, Point _p2)
		: p1(_p1), p2(_p2) {}

	Edge(const Edge& e) {
		p1 = e.p1;
		p2 = e.p2;
	}

	bool operator==(const Edge& e) {
		return (p1 == e.p1 && p2 == e.p2) || 
			(p1 == e.p2 && p2 == e.p1);
	}

	Point p1, p2;
};

/// @brief Status of a triangle - for Triangulation
/// INCOMPLETE: Not yet checked for points overlap
/// COMPLETE: Done checking for points overlap
enum class TriangleStatus {
	INCOMPLETE,
	COMPLETE
};

/// @brief Made of three edges - 
/// Stores the indices of these edges.
/// Indices are of the m_Edges vector in DungeonGenerator
class Triangle {
public:
	Triangle(unsigned a, unsigned b, unsigned c) : 
		 edgeIndexA(a), edgeIndexB(b), edgeIndexC(c),
		 circumCenter(), circumRadiusSquared(0.0f),
		 status(TriangleStatus::INCOMPLETE) {}

	TriangleStatus status;
	unsigned edgeIndexA, edgeIndexB, edgeIndexC;
	Point circumCenter;
	float circumRadiusSquared;
};