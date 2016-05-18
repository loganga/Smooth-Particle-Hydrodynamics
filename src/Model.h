#ifndef __MODEL
#define __MODEL

#define GRAV_ACCEL 980

#include <vector>
#include "glm/glm.hpp"
using namespace std; //makes using vectors easy

class Model
{
public:
	
	Model()
	{
		minBound = glm::vec3(-2.0);
		maxBound = glm::vec3(2.0);
		
		for(int i=0; i<positions.size(); i++)
			positions[i] = positions[i] * 0.5f;
		
	}

	void generatePoints(int amount)
	{
		for (int i = 0; i < amount; i++)
		{
			GLfloat x, y, z;
			x = minBound[0] + static_cast <float> (rand()) / (static_cast <float> (RAND_MAX / (maxBound[0] - minBound[0])));
//			x = std::fmodf((float)rand(), (float)((maxBound[0] - minBound[0] + 1) + minBound[0]));
			y = 0.0f;
//			z = std::fmodf((float)rand(), (float)((maxBound[2] - minBound[2] + 1) + minBound[2]));
			z = minBound[2] + static_cast <float> (rand()) / (static_cast <float> (RAND_MAX / (maxBound[2] - minBound[2])));

			positions.push_back(x * 0.5f);
			positions.push_back(y * 0.5f);
			positions.push_back(z * 0.5f);
			velocity.push_back(0.0f);
			velocity.push_back(0.0f);
			velocity.push_back(0.0f);
			colors.push_back(1.0f);
			colors.push_back(0.0f);
			colors.push_back(1.0f);
			particleCount++;
		}
	}

	void updateParticles(float &deltaT)
	{
		for (int i = 0; i < particleCount; i++) {
			GLfloat x = positions[i * 3];
			GLfloat y = positions[i * 3 + 1];
			GLfloat z = positions[i * 3 + 2];
			GLfloat xVel = velocity[i * 3];
			GLfloat yVel = velocity[i * 3 + 1];
			GLfloat zVel = velocity[i * 3 + 2];
			calculateNewPositions(deltaT, x, y, z, xVel, yVel, zVel);
			calcBoundPositions(x, y, z, xVel, yVel, zVel);
			positions[i * 3] = x;
			positions[i * 3 + 1] = y;
			positions[i * 3 + 2] = z;
			velocity[i * 3] = xVel;
			velocity[i * 3 + 1] = yVel;
			velocity[i * 3 + 2] = zVel;
		}

	}
	
	GLfloat const * getPosition() const
	{ return &positions[0]; }
	
	GLfloat const * getColor() const
	{ return &colors[0]; }
	
	GLuint const * getElements() const
	{ return &elements[0]; }
	
	size_t getVertexCount() const
	{ return positions.size()/3; }
	
	size_t getPositionBytes() const
	{ return positions.size()*sizeof(GLfloat); }
	
	size_t getColorBytes() const
	{ return colors.size()*sizeof(GLfloat); }
	
	size_t getElementBytes() const
	{ return elements.size()*sizeof(GLuint); }

	size_t getParticleCount() const
	{ return particleCount; }
	
	glm::vec3 getCentroid()
	{
		glm::vec3 center = glm::vec3(0);
		float positionCount = 1.0f/(positions.size()/3.0f);
		
		for(int i=0; i<positions.size(); i+=3)
		{
			center[0] += positions[i] * positionCount;
			center[1] += positions[i+1] * positionCount;
			center[2] += positions[i+2] * positionCount;
		}
		
		return center;
	}
	
private:
	glm::vec3 minBound;
	glm::vec3 maxBound;
	vector<GLfloat> positions;
	vector<GLfloat> velocity;
	vector<GLfloat> colors;
	vector<GLuint> elements;
	size_t objectCount;
	size_t particleCount = 0;

	void calculateNewPositions(float &deltaT, GLfloat &x, GLfloat &y, GLfloat &z, GLfloat &xVel, GLfloat &yVel, GLfloat &zVel)
	{
		calcVelocity(deltaT, yVel);
		calcDisplacement(deltaT, y, yVel);
	}

	void calcBoundPositions(GLfloat &x, GLfloat &y, GLfloat &z, GLfloat &xVel, GLfloat &yVel, GLfloat &zVel)
	{
		if (x < minBound[0]) {
			x = minBound[0]; 
			xVel = 0;
		}
		if (x > maxBound[0]) {
			x = maxBound[0];
			xVel = 0;
		}
		if (y < minBound[1]) {
			y = minBound[1]; 
			yVel = 0;
		}
		if (y > maxBound[1]) {
			y = maxBound[1];
			yVel = 0;
		}
		if (z < minBound[2]) {
			z = minBound[2];
			zVel = 0;
		}
		if (z > maxBound[2]) {
			z = maxBound[2];
			zVel = 0;
		}
	}

	void calcVelocity(float &deltaT, GLfloat &vel)
	{
		GLfloat tSquared = deltaT * deltaT;
		vel = vel + (0.5f * -GRAV_ACCEL * tSquared);
	}

	void calcDisplacement(float &deltaT, GLfloat &x, GLfloat &vel)
	{
		x = x + (vel * deltaT);
	}
};

#endif