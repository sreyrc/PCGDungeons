#pragma once

// TODO: Add more to the class if needed
class Point {
public:
	Point() : x(0.0f), y(0.0f), id(INT_MAX) {}

	Point(float _x, float _y, int _id = INT_MAX, bool _isSupra = false)
		: x(_x), y(_y), id(_id), isSuper(_isSupra) {}

	Point(const Point& p) {
		x = p.x;
		y = p.y;
		id = p.id;
		isSuper = p.isSuper;
	}

	float SquaredDistanceTo(const Point& p) const {
		return ((x - p.x) * (x - p.x)) +
			((y - p.y) * (y - p.y));
	}

	Point& operator=(const Point& p) {
		x = p.x;
		y = p.y;
		id = p.id;
		isSuper = p.isSuper;
		return (*this);
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
	unsigned id;
	bool isSuper = false;
};

enum class EdgeStatus {
	AVAILABLE,
	INVALID,
	IN_MST,
	EXTRA
};

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

enum class TriangleStatus {
	INCOMPLETE,
	COMPLETE
};


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