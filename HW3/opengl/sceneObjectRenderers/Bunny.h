#pragma once
#include "Object.h"

class Bunny : public Object {

public:
	float angleZ;
	float angleY;
	float hop;
	float bunnyHopSpeed;
	float hopDirection;
	bool move;
	float x_location;
	float bunnyHoppingRotation;
	bool pressedA ;
	bool pressedD ;


	Bunny();

	void setMatrices(double deltaTime) override;
	void setRotationMatrix(float deltaTime);
	void clearRotationMatrix();
	void setPressedA(bool pressedOrNot);
	void setPressedD(bool pressedOrNot);
    void startGame() override;
};