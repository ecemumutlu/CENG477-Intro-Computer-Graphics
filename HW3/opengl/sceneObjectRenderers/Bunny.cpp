#include "Bunny.h"


#define BUNNY_INITIAL_X 0.0f
#define BUNNY_INITIAL_Y -5.5f
#define BUNNY_INITIAL_HOP_SPEED 25.0f
#define BUNNY_INITIAL_MOVE_SPEED 15.0f
#define BUNNY_INITIAL_ANGLE_Y -90.0f
#define BUNNY_INITIAL_ANGLE_Z 0.0f
#define BUNNY_MAX_HOP_Y (-8.0f)
#define BUNNY_MIN_HOP_Y (-5.5f)

Bunny::Bunny()
	: Object(), angleZ(BUNNY_INITIAL_ANGLE_Z),angleY(BUNNY_INITIAL_ANGLE_Y), hop(BUNNY_INITIAL_Y), bunnyHopSpeed(BUNNY_INITIAL_HOP_SPEED),  
	hopDirection(-1.0f), move(false), bunnyHoppingRotation(M_PI / 2), x_location(BUNNY_INITIAL_X), pressedA(false), pressedD(false)
{
	this->moveSpeed = BUNNY_INITIAL_MOVE_SPEED;
	this->scalingMatrix = glm::scale(glm::mat4(1.0), glm::vec3(1.2, 1.2, 1.2));
	this->rotationMatrix = glm::rotate<float>(glm::mat4(1.0), (-90 / 180.) * M_PI, glm::vec3(0.0, 1.0, 0.0));
}

void Bunny::setMatrices(double deltaTime)  {

	if(this->pressedA){
		float tmp_xloc = this->x_location - (float)  moveSpeed * deltaTime;
		if(tmp_xloc >= -9.5f) this->x_location = tmp_xloc; 
	}
	else if(this->pressedD){
		float tmp_xloc = this->x_location + (float)  moveSpeed * deltaTime;
		if(tmp_xloc <= 9.5f) this->x_location = tmp_xloc;
	}

	// Compute the modeling matrix bunny
	Object::bunnyPosition = glm::vec3(x_location, hop, BUNNY_INITIAL_Z);
	this->translationMatrix = glm::translate(glm::mat4(1.0), Object::bunnyPosition);

	this->setRotationMatrix(deltaTime);

	this->modelingMatrix = this->translationMatrix * this->rotationMatrix * this->scalingMatrix; // starting from right side, rotate around Y to turn back, then rotate around Z some more at each frame, then translate.

	this->bunnyHopSpeed += 0.03 * deltaTime;
	this-> moveSpeed += (0.3) * deltaTime;
	float tmp_hop = hop + hopDirection * bunnyHopSpeed * deltaTime;
	if (tmp_hop < BUNNY_MAX_HOP_Y || tmp_hop > BUNNY_MIN_HOP_Y ) {
		hopDirection *= -1;
	}
	else
		hop = tmp_hop;
}


void Bunny::setRotationMatrix(float deltaTime) { 
	if(this->bunnyTurn == BUNNY_STATIC){
		return;
	}
	else if(this->bunnyTurn == BUNNY_HAPPY){
		float angleRad = (float)(angleY / 180.0) * M_PI;
		this->rotationMatrix = glm::rotate<float>(glm::mat4(1.0), angleRad, glm::vec3(0.0, 1.0, 0.0));
		
		float tmp_angleY = angleY + bunnyHopSpeed ;
		if(angleY >= 270.0f){
			this->clearRotationMatrix();
			angleY = -90.0f;
		}
		else angleY = tmp_angleY;
	}
	else if (this->bunnyTurn == BUNNY_DEAD) {
		float angleRad = (float)(angleZ / 180.0) * M_PI;
		
		this->rotationMatrix = glm::rotate<float>(glm::mat4(1.0), (-90 / 180.) * M_PI, glm::vec3(0.0, 1.0, 0.0));
		this->rotationMatrix *= glm::rotate<float>(glm::mat4(1.0), angleRad, glm::vec3(1.0, 0.0, 0.0));
		//this->rotationMatrix = glm::rotate<float>(glm::mat4(1.0), angleRad, glm::vec3(0.0, 0.0, 1.0));

		float tmp_angleZ= angleZ - bunnyHopSpeed;
		if(angleZ <= -80.0f){
			Object::gameOver = true;
		}
		else angleZ = tmp_angleZ;
	}

}



void Bunny::clearRotationMatrix(){
	this->bunnyTurn = BUNNY_STATIC;
	this->rotationMatrix = glm::rotate<float>(glm::mat4(1.0), (-90 / 180.) * M_PI, glm::vec3(0.0, 1.0, 0.0));
}


void Bunny::setPressedA(bool pressedOrNot){
	this->pressedA = pressedOrNot;
}

void Bunny::setPressedD(bool pressedOrNot){
	this->pressedD = pressedOrNot;
}

void Bunny::startGame() {
	this->findRandomYellowCubeID();
	this->bunnyTurn = BUNNY_STATIC;
	this->rotationMatrix = glm::rotate<float>(glm::mat4(1.0), (-90 / 180.) * M_PI, glm::vec3(0.0, 1.0, 0.0));
	this->moveSpeed = BUNNY_INITIAL_MOVE_SPEED;
	this->angleZ = BUNNY_INITIAL_ANGLE_Z;
	this->angleY = BUNNY_INITIAL_ANGLE_Y;
	this->hop = BUNNY_INITIAL_Y; 
	this->bunnyHopSpeed = BUNNY_INITIAL_HOP_SPEED;
	this->hopDirection = -1.0f;
	this->move = false; 
	this->bunnyHoppingRotation = M_PI / 2; 
	this->x_location = BUNNY_INITIAL_X;
	this->pressedA = false; 
	this->pressedD = false;

}