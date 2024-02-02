#pragma once
#include "Object.h"

class Quad : public Object {
public:
	GLint offsetLoc;
	GLfloat offset;

	Quad(GLfloat offset);
	void setMatrices(double deltaTime)  override;
	void drawObject() override;
	void getUniforms() override;
	void setOffsetLoc(GLint offsetLoc) { this->offsetLoc = offsetLoc;}
    void startGame() override;

};