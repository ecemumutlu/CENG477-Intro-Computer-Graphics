#include "Cube.h"

#define CUBE_RED (glm::vec3(1.0f, 0.0f, 0.0f))
#define CUBE_YELLOW (glm::vec3(1.0f, 1.0f, 0.0f))

#define CUBE_INITIAL_Z -120.0f
#define CUBE_UNVISIBLE_Z (-3.5f)

float Cube::cubeZ = CUBE_INITIAL_Z;
bool Cube::alreadyChanged = false;


Cube::Cube(glm::vec3 translation, int cubeID) : Object(), translationVec(translation), cubeID(cubeID), cubeVanished(false){
	this->moveSpeed = CUBE_QUAD_MOVE_SPEED;
	this->scalingMatrix = glm::scale(glm::mat4(1.0), glm::vec3(1.5f, 3.0f, 1.0f));
	if(this->cubeID == Object::randomYellowCubeID) 
		this->cubeColor = CUBE_YELLOW;
	else
		this->cubeColor = CUBE_RED;
}


void Cube::setMatrices(double deltaTime) {
	
	if(!alreadyChanged){
		if(Cube::cubeZ >= CUBE_UNVISIBLE_Z){
			Cube::cubeZ = CUBE_INITIAL_Z;
		}
		else{
			Cube::cubeZ += moveSpeed * deltaTime;
		}
		
		moveSpeed += CUBE_QUAD_ACCELARATION * deltaTime;
		alreadyChanged = true;
	}
	

	//glm::vec3(-5.0f, -1.0f, -10.f)
	if(this->cubeVanished == false )
		this->translationVec.z = Cube::cubeZ;

	this->translationMatrix = glm::translate(glm::mat4(1.0), this->translationVec);
	this->modelingMatrix = this->translationMatrix  * this->scalingMatrix; // starting from right side, rotate around Y to turn back, then rotate around Z some more at each frame, then translate.

	if(Cube::cubeZ >= -11.5f && this->cubeVanished != true) {
		if(Object::bunnyPosition.x  <= (this->translationVec.x + 2.0f) && Object::bunnyPosition.x  >= (this->translationVec.x - 2.0f)) {
			if(this->cubeColor == CUBE_RED){
				Object::bunnyTurn = BUNNY_DEAD;
				Object::redHit = true;
			}
			else { // CUBE_COLOR = yellow
				Object::bunnyTurn = BUNNY_HAPPY;
				Object::yellowHit = true;
			}
			this->cubeVanished = true;
			this->translationVec.z = 10.0f;
		}
	}

	if(Cube::cubeZ >= CUBE_UNVISIBLE_Z){
		this->translationVec.z = CUBE_INITIAL_Z;
		// reset
		//this->translationVec.z = CUBE_INITIAL_Z;
		this->cubeVanished = false;
		if(this->cubeID == Object::randomYellowCubeID) 
			this->cubeColor = CUBE_YELLOW;
		else
			this->cubeColor = CUBE_RED;
	}

	
}

void Cube::getUniforms() {
	this->modelingMatrixLoc = glGetUniformLocation(this->program.gProgram, "modelingMatrix");
	this->viewingMatrixLoc = glGetUniformLocation(this->program.gProgram, "viewingMatrix");
	this->projectionMatrixLoc = glGetUniformLocation(this->program.gProgram, "projectionMatrix");
	this->eyePosLoc = glGetUniformLocation(this->program.gProgram, "eyePos");
	this->cubeColorPos = glGetUniformLocation(this->program.gProgram, "cubeColor");
}

void Cube::drawObject() {
	glUseProgram(this->program.gProgram);

	glUniformMatrix4fv(this->projectionMatrixLoc, 1, GL_FALSE, glm::value_ptr(this->projectionMatrix));
	glUniformMatrix4fv(this->viewingMatrixLoc, 1, GL_FALSE, glm::value_ptr(this->viewingMatrix));
	glUniformMatrix4fv(this->modelingMatrixLoc, 1, GL_FALSE, glm::value_ptr(this->modelingMatrix));
	glUniform3fv(this->eyePosLoc, 1, glm::value_ptr(this->eyePos));
	glUniform3fv(this->cubeColorPos, 1, glm::value_ptr(this->cubeColor));

	glBindVertexArray(this->vao);
	glBindBuffer(GL_ARRAY_BUFFER, this->gVertexAttribBuffer);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, this->gIndexBuffer);
	glDrawElements(GL_TRIANGLES, this->getIndexDataSize() * 3, GL_UNSIGNED_INT, 0);
	assert(glGetError() == GL_NO_ERROR);

	glBindBuffer(GL_ARRAY_BUFFER, 0); //unbind
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0); //unbind

}

void Cube::startGame() {
	this->translationVec.z = CUBE_INITIAL_Z;
	this->moveSpeed = CUBE_QUAD_MOVE_SPEED;
	Object::redHit = false;
	Object::yellowHit = false;
	Cube::cubeZ = CUBE_INITIAL_Z;
	this->cubeVanished = false;
	if(this->cubeID == Object::randomYellowCubeID) 
		this->cubeColor = CUBE_YELLOW;
	else
		this->cubeColor = CUBE_RED;
}