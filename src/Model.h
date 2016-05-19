#ifndef __MODEL
#define __MODEL

#define GRAV_ACCEL 9.81
#define REF_DENSITY 1000
#define SPEED_OF_SOUND 340.29

#include <vector>
#include "glm/glm.hpp"
using namespace std; //makes using vectors easy

class Model
{
public:
	
	Model()
	{
		minBound = glm::vec3(-1.0);
		maxBound = glm::vec3(1.0);
		
		for(int i=0; i<positions.size(); i++)
			positions[i] = positions[i] * 0.5f;
		
	}

	void generateParticles(int amount)
	{
		for (int i = 0; i < amount; i++)
		{
			GLfloat x, y, z;
			x = minBound[0] + static_cast <float> (rand()) / (static_cast <float> (RAND_MAX / (maxBound[0] - minBound[0])));
			y = minBound[1] + static_cast <float> (rand()) / (static_cast <float> (RAND_MAX / (maxBound[1] - minBound[1])));
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
			acceleration.push_back(0.0f);
			acceleration.push_back(-GRAV_ACCEL);
			acceleration.push_back(0.0f);
			density.push_back(0.0f);
			mass = 1.0f;
			radius = 1.0f;
			viscosity = 0.1f;
//			bulkModulus = (SPEED_OF_SOUND * SPEED_OF_SOUND) * REF_DENSITY;
			bulkModulus = 1000;
			particleCount++;

		}
//		calcDensity();
	}

	void updateParticles(float &deltaT, glm::mat4 rot)
	{
		for (int i = 0; i < particleCount; i++) {
			GLfloat x = positions[i * 3];
			GLfloat y = positions[i * 3 + 1];
			GLfloat z = positions[i * 3 + 2];
			GLfloat xVel = velocity[i * 3];
			GLfloat yVel = velocity[i * 3 + 1];
			GLfloat zVel = velocity[i * 3 + 2];
			calcDensity();
			calcAcceleration(deltaT);
			glm::vec4 temp = glm::vec4(acceleration[i * 3], acceleration[i * 3 + 1], acceleration[i * 3 + 2], 0) * rot;
			acceleration[i * 3] = temp[0];
			acceleration[i * 3 + 1] = temp[1];
			acceleration[i * 3 + 2] = temp[2];
			calculateNewPositions(deltaT, i, x, y, z, xVel, yVel, zVel);
			positions[i * 3] = x;
			positions[i * 3 + 1] = y;
			positions[i * 3 + 2] = z;
			calcBoundPositions(i, x, y, z);
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
	vector<GLfloat> acceleration;
	vector<GLfloat> colors;
	vector<GLuint> elements;
	vector<GLfloat> density;
	GLfloat viscosity;
	GLfloat radius;
	GLfloat mass;
	GLfloat bulkModulus;
	size_t objectCount;
	size_t particleCount = 0;

	void calculateNewPositions(float &deltaT, int i, GLfloat &x, GLfloat &y, GLfloat &z, GLfloat &xVel, GLfloat &yVel, GLfloat &zVel)
	{
		calcVelocity(deltaT, i);
		calcDisplacement(deltaT, x, xVel);
		calcDisplacement(deltaT, y, yVel);
		calcDisplacement(deltaT, z, zVel);
	}

	void calcBoundPositions(int particleNum, GLfloat &x, GLfloat &y, GLfloat &z)
	{
		if (x < minBound[0]) {
			reflect(particleNum, 0, minBound[0]);
		}
		if (x > maxBound[0]) {
			reflect(particleNum, 0, maxBound[0]);
		}
		if (y < minBound[1]) {
			reflect(particleNum, 1, minBound[1]);
		}
		if (y > maxBound[1]) {
			reflect(particleNum, 1, maxBound[1]);
		}
		if (z < minBound[2]) {
			reflect(particleNum, 2, minBound[2]);
		}
		if (z > maxBound[2]) {
			reflect(particleNum, 2, maxBound[2]);
		}
	}

	void calcVelocity(float &deltaT, int i)
	{
		GLfloat tSquared = deltaT * deltaT;
		velocity[i * 3] += (0.5f * acceleration[i * 3] * tSquared);
		velocity[i * 3 + 1] += (0.5f * acceleration[i * 3 + 1] * tSquared);
		velocity[i * 3 + 2] += (0.5f * acceleration[i * 3 + 2] * tSquared);
	}

	void calcDisplacement(float &deltaT, GLfloat &x, GLfloat &vel)
	{
		x = x + (vel * deltaT);
	}

	void calcDensity()
	{
		GLfloat coefficient = (4 * mass) / (PI * powf(radius, 8.0f));
		GLfloat sum  = 0;
		for (int i = 0; i < particleCount - 1; i++)
		{
			density[i] += (4 * mass) / (PI * (radius * radius));
			for (int j = i + 1; j < particleCount; j++)
			{
				GLfloat dx	= positions[3 * i] - positions[3 * j];
				GLfloat dy	= positions[3 * i + 1] - positions[3 * j + 1];
				GLfloat dz	= positions[3 * i + 2] - positions[3 * j + 2];
				GLfloat r2	= dx*dx + dy*dy + dz*dz;
				GLfloat z = (radius * radius) - r2;
				if (z > 0)
				{
					GLfloat rho = coefficient*z*z*z;
					density[i] += rho;
					density[j] += rho;
				}
			}
		}
	}

	void calcAcceleration(GLfloat deltaT)
	{
		for (int i = 0; i < particleCount; i++)
		{
			acceleration[3 * i] = 0;
			acceleration[3 * i + 1] = -GRAV_ACCEL / deltaT;
			acceleration[3 * i + 2] = 0;
		}

		GLfloat C0 = mass / PI / powf(radius, 4.0f);
		GLfloat Cp = 15 * bulkModulus;
		GLfloat Cv = -40 * viscosity;

		for (int i = 0; i < particleCount - 1; i++)
		{
			GLfloat dens_i = density[i];
			for (int j = i + 1; j < particleCount; j++)
			{
				GLfloat dx = positions[3 * i];
				GLfloat dy = positions[3 * i + 1];
				GLfloat dz = positions[3 * i + 2];
				GLfloat r2 = dx * dx + dy * dy + dz * dz;
				if (r2 < (radius * radius))
				{
					GLfloat dens_j = density[j];
					GLfloat q = sqrt(r2) / radius;
					GLfloat u = 1 - q;
					GLfloat w0 = C0 * u / dens_i / dens_j;
					GLfloat wp = w0 * Cp * (dens_i + dens_j - 2 * REF_DENSITY) * u/q;
					GLfloat wv = w0 * Cv;
					GLfloat dvx = velocity[3 * i] - velocity[3 * j];
					GLfloat dvy = velocity[3 * i + 1] - velocity[3 * j + 1];
					GLfloat dvz = velocity[3 * i + 2] - velocity[3 * j + 2];
					acceleration[3 * i] += wp * dx + wv * dvx;
					acceleration[3 * i + 1] += wp * dy + wv * dvy;
					acceleration[3 * i + 2] += wp * dz + wv * dvz;
					acceleration[3 * j] -= wp * dx + wv * dvx;
					acceleration[3 * j + 1] -= wp * dy + wv * dvy;
					acceleration[3 * j + 2] -= wp * dz + wv * dvz;
				}
			}
		}
	}

	void reflect(int particleNum, int dimension, GLfloat bound)
	{
		// dampening constant
		const GLfloat dampen = 0.5f;
		int i = particleNum * 3 + dimension;

		//Scale the distance traveled based on collision time
		if (velocity[i] == 0)
			return;

		GLfloat tbounce = (positions[i] - bound) / velocity[i];
		positions[particleNum * 3] -= velocity[particleNum * 3] * (1 - dampen) * tbounce;
		positions[particleNum * 3 + 1] -= velocity[particleNum * 3 + 1] * (1 - dampen) * tbounce;
		positions[particleNum * 3 + 2] -= velocity[particleNum * 3 + 2] * (1 - dampen) * tbounce;

		positions[i] = 2 * bound - positions[i];
		velocity[i] = -velocity[i];

		velocity[particleNum * 3] *= dampen;
		velocity[particleNum * 3 + 1] *= dampen;
		velocity[particleNum * 3 + 2] *= dampen;
	}

};

#endif