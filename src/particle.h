// particle.h

#include "includes.h"

// struct Module
// {
// 	virtual double norm();
// 	virtual Module operator+(const Module& other);
// 	virtual Module operator-(const Module& other);
// 	virtual Module operator*=(double scale);
// 	virtual Module operator*(const double& scale, Module rhs)
// 	{
// 		return rhs *= scale;
// 	}

// };

// struct Point : Module
// {
// public:
// 	double x;

// 	Point(double x = 0)
// 	{
// 		this->x = x;
// 	}

// 	double norm()
// 	{
// 		return abs(x);
// 	}

// 	Point operator+(const Point& other)
// 	{
// 		return Point(this->x + other.x);
// 	}

// 	Point operator-(const Point& other)
// 	{
// 		return Point(this->x - other.x);
// 	}

// 	Point operator*=(double scale)
// 	{
// 		this->x *= scale;
// 	}
// };

// struct Pair : Module
// {
// public:
// 	double x, y;

// 	Pair(double x = 0, double y = 0)
// 	{
// 		this->x = x;
// 		this->y = y;
// 	}

// 	double norm()
// 	{
// 		return sqrt(x * x + y * y);
// 	}

// 	Pair operator+(const Pair& other)
// 	{
// 		return Pair(this->x + other.x, this->y + other.y);
// 	}

// 	Pair operator-(const Pair& other)
// 	{
// 		return Pair(this->x - other.x, this->y - other.y);
// 	}

// 	Pair operator*=(double scale)
// 	{
// 		this->x *= scale;
// 		this->y *= scale;
// 	}
// };

// struct Triple : Module
// {
// public:
// 	double x, y, z;

// 	Triple(double x = 0, double y = 0, double z = 0)
// 	{
// 		this->x = x;
// 		this->y = y;
// 		this->z = z;
// 	}

// 	double norm()
// 	{
// 		return sqrt(x * x + y * y + z * z);
// 	}

// 	Triple operator+(const Triple& other)
// 	{
// 		return Triple(this->x + other.x, this->y + other.y, this->z + other.z);
// 	}

// 	Triple operator-(const Triple& other)
// 	{
// 		return Triple(this->x - other.x, this->y - other.y, this->z - other.z);
// 	}

// 	Triple operator*=(double scale)
// 	{
// 		this->x *= scale;
// 		this->y *= scale;
// 		this->z *= scale;
// 	}
// };



// struct Point 
// {
// public:
// 	double x;

// 	Point(double x = 0)
// 	{
// 		this->x = x;
// 	}

// 	double norm()
// 	{
// 		return abs(x);
// 	}

// 	Point operator+(const Point& other)
// 	{
// 		return Point(this->x + other.x);
// 	}

// 	Point operator-(const Point& other)
// 	{
// 		return Point(this->x - other.x);
// 	}
// };

// Point operator*(const double& scale, Point rhs)
// {
// 	return Point(rhs.x * scale);
// }

// Point operator/(const double& scale, Point rhs)
// {
// 	return Point(rhs.x / scale);
// }

// Point operator*(Point lhs, const double& scale)
// {
// 	return scale * lhs;
// }

// Point operator/(Point lhs, const double& scale)
// {
// 	return lhs / scale;
// }


// // position, velocity, and acceleration take values in a double-module M
// template <typename M> class Particle
// {
// public:
// 	double mass;
// 	M position;
// 	M velocity;
// 	M acceleration;
// 	double drag;

// 	Particle(double m, double x, double v = 0, double a = 0, double d = 1)
// 	{
// 		mass = m;
// 		position = x;
// 		velocity = v;
// 		acceleration = a;
// 		drag = relaxation(d);
// 	}

// 	void tick()
// 	{
// 		velocity += acceleration / SR;
// 		velocity *= drag;
// 		position += velocity / SR;
// 	}

// 	void prepare()
// 	{
// 		acceleration = 0;
// 	}

// 	void pull(M force)
// 	{
// 		acceleration += force / mass;
// 	}
// };

// template <typename M> class Spring
// {
// public:
// 	Particle<M>* first;
// 	Particle<M>* second;
// 	double k;
// 	double equilibrium;

// 	Spring(Particle<M>* first, Particle<M>* second, double k, double equilibrium)
// 	{
// 		this->first = first;
// 		this->second = second;
// 		this->k = k;
// 		this->equilibrium = equilibrium;
// 	}

// 	void tick()
// 	{
// 		double distance = (second->position - first->position).norm();
// 		double displacement = distance - equilibrium;
// 		first->prepare();
// 		second->prepare();
// 		first->pull(k * displacement * (second->position - first->position) / distance);
// 		second->pull(-k * displacement * (first->position - second->position) / distance);
// 	}
// };


// position, velocity, and acceleration take values in a double-module M
class Particle
{
public:
	double mass;
	double position;
	double velocity;
	double acceleration;
	double drag;

	void initialize(double m, double x, double v = 0, double a = 0, double d = 1)
	{
		mass = m;
		position = x;
		velocity = v;
		acceleration = a;
		drag = relaxation(d);
	}

	Particle();

	// Particle(double m, double x, double v = 0, double a = 0, double d = 1)
	// {
	// 	mass = m;
	// 	position = x;
	// 	velocity = v;
	// 	acceleration = a;
	// 	drag = relaxation(d);
	// }

	void reset(double x)
	{
		position = x;
		velocity = 0;
		acceleration = 0;
	}

	void tick()
	{
		velocity += acceleration / SR;
		velocity *= drag;
		position += velocity / SR;
	}

	void prepare()
	{
		acceleration = 0;
	}

	void pull(double force)
	{
		acceleration += force / mass;
	}
};

class Spring
{
public:
	Particle* first;
	Particle* second;
	double k;
	double equilibrium;

	void initialize(Particle* first, Particle* second, double k, double equilibrium)
	{
		this->first = first;
		this->second = second;
		this->k = k;
		this->equilibrium = equilibrium;
	}

	Spring();

	// Spring(Particle* first, Particle* second, double k, double equilibrium)
	// {
	// 	this->first = first;
	// 	this->second = second;
	// 	this->k = k;
	// 	this->equilibrium = equilibrium;
	// }

	void tick()
	{
		double distance = second->position - first->position;
		double displacement = distance - equilibrium;
		first->pull(k * displacement);
		second->pull(-k * displacement);
	}
};

class Gravity
{
public:
	Particle* first;
	Particle* second;
	double G;
	double epsilon;

	void initialize(Particle* first, Particle* second, double G, double epsilon = 0.0001)
	{
		this->first = first;
		this->second = second;
		this->G = G;
		this->epsilon = epsilon;
	}

	Gravity();

	// Gravity(Particle* first, Particle* second, double G, double epsilon =  0.0001)
	// {
	// 	this->first = first;
	// 	this->second = second;
	// 	this->G = G;
	// 	this->epsilon = epsilon;
	// }

	void tick()
	{
		double distance = second->position - first->position;
		double cubed = abs(distance * distance * distance);
		first->pull(G * first->mass * second->mass * distance / (epsilon + cubed));
		first->pull(-G * first->mass * second->mass * distance / (epsilon + cubed));
	}
};