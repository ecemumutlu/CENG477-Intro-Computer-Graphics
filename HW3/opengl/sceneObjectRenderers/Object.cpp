#include "Object.h"
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>
#include <ctime>
using namespace std;



void Object::setUniforms(GLint projectionMatrixLoc, GLint viewingMatrixLoc, GLint modelingMatrixLoc, GLint eyePosLoc) {
	glUniformMatrix4fv(projectionMatrixLoc, 1, GL_FALSE, glm::value_ptr(this->projectionMatrix));
	glUniformMatrix4fv(viewingMatrixLoc, 1, GL_FALSE, glm::value_ptr(this->viewingMatrix));
	glUniformMatrix4fv(modelingMatrixLoc, 1, GL_FALSE, glm::value_ptr(this->modelingMatrix));
	glUniform3fv(eyePosLoc, 1, glm::value_ptr(this->eyePos));

	//glm::vec3* color_ins = new glm::vec3(255.0f / 255.0f, 250.0f / 255.0f, 71.0f / 255.0f);
	//glUniform3f(colorLocation, color_ins->r, color_ins->g, color_ins->b);
	//delete color_ins;
}

int Object::randomYellowCubeID = 1;
glm::vec3 Object::bunnyPosition = glm::vec3(0.0f, -4.0f, -7.0f);
bool Object::gameOver = false;
int Object::bunnyTurn = BUNNY_STATIC;
bool Object::yellowHit = false;
bool Object::redHit = false;


void Object::drawObject() {
	glUseProgram(this->program.gProgram);

	glUniformMatrix4fv(this->projectionMatrixLoc, 1, GL_FALSE, glm::value_ptr(this->projectionMatrix));
	glUniformMatrix4fv(this->viewingMatrixLoc, 1, GL_FALSE, glm::value_ptr(this->viewingMatrix));
	glUniformMatrix4fv(this->modelingMatrixLoc, 1, GL_FALSE, glm::value_ptr(this->modelingMatrix));
	glUniform3fv(this->eyePosLoc, 1, glm::value_ptr(this->eyePos));

	glBindVertexArray(this->vao);
	glBindBuffer(GL_ARRAY_BUFFER, this->gVertexAttribBuffer);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, this->gIndexBuffer);
	glDrawElements(GL_TRIANGLES, this->getIndexDataSize() * 3, GL_UNSIGNED_INT, 0);
	assert(glGetError() == GL_NO_ERROR);

	glBindBuffer(GL_ARRAY_BUFFER, 0); //unbind
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0); //unbind

}

void Object::findRandomYellowCubeID(){
	std::srand(static_cast<unsigned int>(std::time(nullptr)));

    // Generate a random value between 0 and 2
	Object::randomYellowCubeID = std::rand() % 3;
    return ;
}



int Object::getVertexDataSize(bool inBytes) { 
	if (!inBytes) 
		return oVertices.size(); 
	else  
		return oVertices.size() * 3 * sizeof(GLfloat); 
}

int Object::getNormalDataSize(bool inBytes) {
	if (!inBytes) 
		return oNormals.size();
	else 
		return oNormals.size() * 3 * sizeof(GLfloat);
}

int Object::getColorDataSize(bool inBytes) {
	if (!inBytes) 
		return oColors.size();
	else 
		return oColors.size() * 3 * sizeof(GLfloat);
}

int Object::getIndexDataSize(bool inBytes) {
	if (!inBytes) 
		return oFaces.size();
	else 
		return oFaces.size() * 3 * sizeof(GLuint);
}

int Object::getTextureDataSize(bool inBytes){
	if (!inBytes) 
		return oTextures.size();
	else 
		return oTextures.size() * 2 * sizeof(GLfloat);
}


void Object::getUniforms() {
	this->modelingMatrixLoc = glGetUniformLocation(this->program.gProgram, "modelingMatrix");
	this->viewingMatrixLoc = glGetUniformLocation(this->program.gProgram, "viewingMatrix");
	this->projectionMatrixLoc = glGetUniformLocation(this->program.gProgram, "projectionMatrix");
	this->eyePosLoc = glGetUniformLocation(this->program.gProgram, "eyePos");
}










