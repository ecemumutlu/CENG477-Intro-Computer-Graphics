#include "Quad.h"

Quad::Quad(float offset) 
	: Object(), offset(offset)
{
	this->moveSpeed = CUBE_QUAD_MOVE_SPEED;
	this->rotationMatrix = glm::rotate<float>(glm::mat4(1.0), (-90 / 180.) * M_PI, glm::vec3(1.0, 0.0, 0.0));
	this->scalingMatrix = glm::scale(glm::mat4(1.0), glm::vec3(17.0f, 1000.0f, 1.0f)); //because we rotate around x, the alongated look towards z-axis is achieved by applying scale more on y-axis

}

void Quad::setMatrices(double deltaTime)  {
	//glm::vec3(-5.0f, -1.0f, -10.f)
	this->translationMatrix = glm::translate(glm::mat4(1.0), glm::vec3(0.0f, -10.0f, 0.0f));
	this->modelingMatrix = this->translationMatrix * this->rotationMatrix * this->scalingMatrix;  // starting from right side, rotate around Y to turn back, then rotate around Z some more at each frame, then translate.
	this->offset -= moveSpeed * deltaTime;

	moveSpeed += CUBE_QUAD_ACCELARATION * deltaTime;
}


void Quad::drawObject() {
	glUseProgram(this->program.gProgram);
	glUniformMatrix4fv(this->projectionMatrixLoc, 1, GL_FALSE, glm::value_ptr(this->projectionMatrix));
	glUniformMatrix4fv(this->viewingMatrixLoc, 1, GL_FALSE, glm::value_ptr(this->viewingMatrix));
	glUniformMatrix4fv(this->modelingMatrixLoc, 1, GL_FALSE, glm::value_ptr(this->modelingMatrix));
	glUniform3fv(this->eyePosLoc, 1, glm::value_ptr(this->eyePos));
	glUniform1f(this->offsetLoc, this->offset);

	glBindVertexArray(this->vao);
	glBindBuffer(GL_ARRAY_BUFFER, this->gVertexAttribBuffer);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, this->gIndexBuffer);
	glDrawElements(GL_TRIANGLES, this->getIndexDataSize() * 3, GL_UNSIGNED_INT, 0);

	glBindBuffer(GL_ARRAY_BUFFER, 0); //unbind
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0); //unbind

}


void Quad::getUniforms(){

	this->modelingMatrixLoc = glGetUniformLocation(this->program.gProgram, "modelingMatrix");
	this->viewingMatrixLoc = glGetUniformLocation(this->program.gProgram, "viewingMatrix");
	this->projectionMatrixLoc = glGetUniformLocation(this->program.gProgram, "projectionMatrix");
	this->eyePosLoc = glGetUniformLocation(this->program.gProgram, "eyePos");
	this->offsetLoc = glGetUniformLocation(this->program.gProgram, "offset");
}

void Quad::startGame() {
	this->moveSpeed = CUBE_QUAD_MOVE_SPEED;
	this->offset = 70.0f;
}