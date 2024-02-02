#pragma once
#include "Object.h"

class Cube : public Object {
public:
	static float cubeZ;
	static bool alreadyChanged ;

	bool cubeVanished;
	glm::vec3 translationVec;
	float z_location;
	GLint cubeColorPos;
	int cubeID;
	glm::vec3 cubeColor;
	bool vanishCube;


	Cube(glm::vec3 translation, int cubeID);

	void setMatrices(double deltaTime) override;
	void getUniforms() override;
	void drawObject() override;
    void startGame() override;
};