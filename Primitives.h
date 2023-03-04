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
	Triangle(Point pA, Point pB, Point pC) : status(Status::INCOMPLETE) {

		vA = pA;
		float temp = (pB.x - pA.x) * (pC.y - pA.y) - (pC.x - pA.x) * (pB.y - pA.y);
		bool isCounterClockwise = temp > 0;

		vA = pA;
		vB = isCounterClockwise ? pB : pC;
		vC = isCounterClockwise ? pC : pB;

		circumCenter = CalculateCircumCenter();
		circumRadiusSquared = vA.SquaredDistanceTo(circumCenter);
	}

	Point CalculateCircumCenter() {

		// Formulae: https://tinyurl.com/yvsw7s5f

		Point SqrA(powf(vA.x, 2), powf(vA.y, 2));
		Point SqrB(powf(vB.x, 2), powf(vB.y, 2));
		Point SqrC(powf(vC.x, 2), powf(vC.y, 2));

		float D = (vA.x * (vB.y - vC.y) + vB.x * (vC.y - vA.y) + vC.x * (vA.y - vB.y)) * 2.0f;
		float x = ((SqrA.x + SqrA.y) * (vB.y - vC.y) + (SqrB.x + SqrB.y) * 
			(vC.y - vA.y) + (SqrC.x + SqrC.y) * (vA.y - vB.y)) / D;
		float y = ((SqrA.x + SqrA.y) * (vC.x - vB.x) + (SqrB.x + SqrB.y) * 
			(vA.x - vC.x) + (SqrC.x + SqrC.y) * (vB.x - vA.x)) / D;

		return Point(x, y);
	}

	Point GetCircumCenter() { return circumCenter; }
	float GetCircumRadiusSquared() { return circumRadiusSquared; }

	Status status;
	Point vA, vB, vC;
private:
	Point circumCenter;
	float circumRadiusSquared;
	// TODO: Store edges as well
};