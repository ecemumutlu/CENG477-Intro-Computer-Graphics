#include "Sky.h"
#include <iostream>
void Sky::getUniforms(){
    //this->bindTexture();

    this->textureLoc = glGetUniformLocation(this->program.gProgram, "u_Texture");
}

void Sky::setMatrices(double deltaTime) {

    
}


void Sky::drawObject() {
	glUseProgram(this->program.gProgram);
	this->bindTexture();
    
    glUniform1i(this->textureLoc, 0);
	
	glBindVertexArray(this->vao);
	glBindBuffer(GL_ARRAY_BUFFER, this->gVertexAttribBuffer);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, this->gIndexBuffer);
	glDrawElements(GL_TRIANGLES, this->getIndexDataSize() * 3, GL_UNSIGNED_INT, 0);

	glBindBuffer(GL_ARRAY_BUFFER, 0); //unbind
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0); //unbind

}



void Sky::bindTexture(){
    this->texture.Bind(this->slot);
}