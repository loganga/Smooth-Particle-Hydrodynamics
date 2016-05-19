#ifndef __PARTICLE
#define __PARTICLE
#include "src/glm/detail/type_vec3.hpp"
#include "src/gl3w/glcorearb.h"

class Particle
{
public:
	glm::vec3 pos, vel, acc;
	GLfloat r, g, b, a, size, density;

	Particle() {}

	Particle(GLfloat x, GLfloat y, GLfloat z, GLfloat r, GLfloat g, GLfloat b, GLfloat a)
	{
		pos[0] = x;
		pos[1] = y;
		pos[2] = z;
		vel = glm::vec3(0.0f);
		this->r = r;
		this->g = g;
		this->b = b;
		this->a = a;
		size = 0.1f;
	}

	void setPos(GLfloat x, GLfloat y, GLfloat z)
	{
		pos[0] = x;
		pos[1] = y;
		pos[2] = z;
	}

};

#endif