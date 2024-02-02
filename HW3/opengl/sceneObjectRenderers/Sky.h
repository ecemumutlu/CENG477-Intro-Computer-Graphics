#pragma once
#include "Object.h"
#include "Texture.h"

class Sky: public Object {
public:
    Sky(const std::string& path, int slot) : Object(), texture(path), slot(slot) {}
    ~Sky();

    GLint textureLoc;
    TextureSky texture;
    int slot;

	void setMatrices(double deltaTime)  override;
	void drawObject() override;
	void getUniforms() override;

    void bindTexture();
};