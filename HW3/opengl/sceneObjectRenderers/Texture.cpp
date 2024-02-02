#include "Texture.h"
#define STB_IMAGE_IMPLEMENTATION
#include "../stb_image/stb_image.h"
#include <iostream>

TextureSky::TextureSky(const std::string& path) 
	: m_RendererID(0), m_FilePath(path), m_LocalBuffer(nullptr), 
	m_Width(0), m_Height(0), m_BPP(0)
{
	stbi_set_flip_vertically_on_load(1);
	m_LocalBuffer = stbi_load(path.c_str(), &m_Width, &m_Height, &m_BPP, 4); //desired channels = 4 means RGBA

	glGenTextures(1, &m_RendererID);
	glBindTexture(GL_TEXTURE_2D, m_RendererID);//bind

	//below 4 are important
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR );
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR );
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);

	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA8, m_Width, m_Height, 0, GL_RGBA, GL_UNSIGNED_BYTE, m_LocalBuffer);
	glBindTexture(GL_TEXTURE_2D, 0); //unbind

	if (m_LocalBuffer)
		stbi_image_free(m_LocalBuffer);
}

TextureSky::~TextureSky() {
	glDeleteTextures(1, &m_RendererID);
}

void TextureSky::Bind(unsigned int slot ) const {
	glActiveTexture(GL_TEXTURE0 ); //in opengl there is a limitation on the number of textures
	glBindTexture(GL_TEXTURE_2D, m_RendererID);//bind
}

void TextureSky::Unbind() const {
	glBindTexture(GL_TEXTURE_2D, 0); //unbind
}