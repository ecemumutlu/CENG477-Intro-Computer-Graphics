
#pragma once
#include "GL/glew.h"
#include <vector>
#include <glm/glm.hpp> // GL Math library header
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp> 
#include <math.h>
#include <cmath>
#include <iostream>
#include "Program.h"

#define M_PI       3.14159265358979323846   // pi
#define BUNNY_TYPE 0
#define CUBE_TYPE 1
#define QUAD_TYPE 2
#define SKY_TYPE 3

//positions of objects
#define MOVE_LEFT -1
#define NO_MOVE 0
#define MOVE_RIGHT 1

#define BUNNY_STATIC 0
#define BUNNY_DEAD 1
#define BUNNY_HAPPY 2

#define CUBE_QUAD_MOVE_SPEED 45.0f
#define CUBE_QUAD_ACCELARATION (0.5f)
#define BUNNY_INITIAL_Z -11.0f

struct Vertex
{
	Vertex(GLfloat inX, GLfloat inY, GLfloat inZ) : x(inX), y(inY), z(inZ) { }
	GLfloat x, y, z;
};


struct Texture
{
	Texture(GLfloat inU, GLfloat inV) : u(inU), v(inV) { }
	GLfloat u, v;
};

struct Normal
{
	Normal(GLfloat inX, GLfloat inY, GLfloat inZ) : x(inX), y(inY), z(inZ) { }
	GLfloat x, y, z;
};

struct Face
{
	Face(int v[], int t[], int n[]) {
		vIndex[0] = v[0];
		vIndex[1] = v[1];
		vIndex[2] = v[2];
		tIndex[0] = t[0];
		tIndex[1] = t[1];
		tIndex[2] = t[2];
		nIndex[0] = n[0];
		nIndex[1] = n[1];
		nIndex[2] = n[2];
	}
	GLuint vIndex[3], tIndex[3], nIndex[3];
};



class Object {
public:
	static int randomYellowCubeID;
	static glm::vec3 bunnyPosition;
	static bool gameOver;
	static int bunnyTurn;
	static bool yellowHit;
	static bool redHit;

	std::vector<Vertex> oVertices;
	std::vector<glm::vec3> oColors;
	std::vector<Texture> oTextures;
	std::vector<Normal> oNormals;
	std::vector<Face> oFaces;
	int object_type;
	glm::mat4 translationMatrix;
	glm::mat4 scalingMatrix;
	glm::mat4 rotationMatrix;
	glm::mat4 modelingMatrix;
	glm::mat4 projectionMatrix;
	glm::mat4 viewingMatrix;
	glm::vec3 eyePos;
	float moveSpeed;
	

	GLuint vao;
	GLuint gVertexAttribBuffer;
	GLuint gIndexBuffer;

	Program program;
	GLint modelingMatrixLoc; GLint viewingMatrixLoc; GLint projectionMatrixLoc; GLint eyePosLoc;

	Object()
		: translationMatrix(glm::mat4(1.0)),
		scalingMatrix(glm::mat4(1.0)),
		rotationMatrix(glm::mat4(1.0)),
		modelingMatrix(glm::mat4(1.0)),
		projectionMatrix(glm::mat4(1.0)),
		viewingMatrix(glm::mat4(1.0)),
		eyePos(glm::vec3(0, 0, 0)),
		vao(0), gVertexAttribBuffer(0), gIndexBuffer(0), program(),
		modelingMatrixLoc(-1), viewingMatrixLoc(-1), projectionMatrixLoc(-1), eyePosLoc(-1),
		moveSpeed(0.0f)
	{
		findRandomYellowCubeID();
	}

	// Constructor
	Object(const std::vector<Vertex>& vertices,
		const std::vector<glm::vec3>& colors,
		const std::vector<Texture>& textures,
		const std::vector<Normal>& normals,
		const std::vector<Face>& faces,
		int obj_type
		)
		: oVertices(vertices), oColors(colors), oTextures(textures), oNormals(normals), oFaces(faces), object_type(obj_type), 
		translationMatrix(glm::mat4(1.0)), 
		scalingMatrix(glm::mat4(1.0)), 
		rotationMatrix(glm::mat4(1.0)), 
		modelingMatrix(glm::mat4(1.0)), 
		projectionMatrix(glm::mat4(1.0)), 
		viewingMatrix(glm::mat4(1.0)), 
		eyePos(glm::vec3(0,0,0)),
		vao(0), gVertexAttribBuffer(0), gIndexBuffer(0), program(),
		modelingMatrixLoc(-1), viewingMatrixLoc(-1), projectionMatrixLoc(-1), eyePosLoc(-1),
		moveSpeed(0.0f)
	{
		findRandomYellowCubeID();
	}

	~Object() {}

	int getVertexDataSize(bool inBytes = false);
	int getNormalDataSize(bool inBytes = false);
	int getColorDataSize(bool inBytes = false);
	int getIndexDataSize(bool inBytes = false);
	int getTextureDataSize(bool inBytes = false);

	inline void setProjectionMatrix(glm::mat4 proj) { this->projectionMatrix = proj; }
	inline void setViewingMatrix(glm::mat4 view) { this->viewingMatrix = view; }
	inline void setObjectType(int obj_type) { this->object_type = obj_type; };

	void setUniforms(GLint projectionMatrixLoc, GLint viewingMatrixLoc, GLint modelingMatrixLoc, GLint eyePosLoc);
	void findRandomYellowCubeID();

	virtual void drawObject();
	virtual void setMatrices(double deltaTime) {}
	virtual void getUniforms();
	virtual void startGame() {};



};

