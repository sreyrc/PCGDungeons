#pragma once

// TODO: Add more to the class if needed
class Point {
public:
	Point() : x(0.0f), y(0.0f) {}

	Point(float _x, float _y)
		: x(_x), y(_y) {}

	Point(const Point& p) {
		x = p.x;
		y = p.y;
		isSupra = p.isSupra;
	}

	float SquaredDistanceTo(const Point& p) const {
		return ((x - p.x) * (x - p.x)) +
			((y - p.y) * (y - p.y));
	}

	Point& operator=(const Point& p) {
		x = p.x;
		y = p.y;
		isSupra = p.isSupra;
		return (*this);
	}

	bool operator==(const Point& p) {
		return (fabs(x - p.x) < FLT_EPSILON)
			&& (fabs(y - p.y) < FLT_EPSILON);
	}

	float x, y;
	bool isSupra = false;
};

// TODO: Add more to the class if needed
class Edge {
public:
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

enum class Status {
	INCOMPLETE,
	COMPLETE
};


class Triangle {
public:
	Triangle(unsigned a, unsigned b, unsigned c) : 
		 edgeIndexA(a), edgeIndexB(b), edgeIndexC(c),
		 circumCenter(), circumRadiusSquared(0.0f),
		 status(Status::INCOMPLETE) {}

	Status status;
	unsigned edgeIndexA, edgeIndexB, edgeIndexC;
	Point circumCenter;
	float circumRadiusSquared;
};